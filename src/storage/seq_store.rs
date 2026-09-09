//! 2-bit encoded genomic sequence storage organized by taxonomy and assembly.
//!
//! # Format Specifications
//!
//! - **v2 Format (`SCIESSEQ`)**:
//!   - 64-byte aligned header with magic bytes `SCIESSEQ`, version, length, checksum,
//!     and ambiguity payload descriptors.
//!   - Packed 2-bit canonical sequence payload (4 bases per byte: `A=00, C=01, G=10, T=11`).
//!   - Ambiguity sidecar table preserving IUPAC degenerate bases and RNA Uracil losslessly.
//!   - CRC32 payload checksum guarding against disk truncation and bit rot.
//!   - Atomic write-and-rename guarantees crash consistency.
//!
//! - **v1 Legacy Format**:
//!   - 8-byte length prefix followed by raw 2-bit payload (supported for read-only backward compatibility).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use memmap2::MmapOptions;

use crate::{BioDbError, Result};

/// Format magic bytes for `.seq v2`.
pub const SEQ_V2_MAGIC: &[u8; 8] = b"SCIESSEQ";

/// Current format version.
pub const SEQ_V2_VERSION: u16 = 2;

/// Fixed 64-byte header size for v2 files.
pub const HEADER_SIZE_V2: usize = 64;

static TEMP_FILE_COUNTER: AtomicU64 = AtomicU64::new(1);

/// Standard IEEE 802.3 CRC32 implementation without external dependencies.
pub fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &byte in data {
        crc ^= byte as u32;
        for _ in 0..8 {
            let mask = -((crc & 1) as i32) as u32;
            crc = (crc >> 1) ^ (0xEDB8_8320 & mask);
        }
    }
    !crc
}

/// Sanitizes path components to prevent directory traversal vulnerabilities.
pub fn sanitize_component<'a>(name: &'a str, context: &str) -> Result<&'a str> {
    if name.is_empty() {
        return Err(BioDbError::PathTraversal(format!(
            "{context} component cannot be empty"
        )));
    }
    if name.contains("..") || name.contains('/') || name.contains('\\') || name.contains('\0') {
        return Err(BioDbError::PathTraversal(format!(
            "Forbidden characters in {context} name '{name}'"
        )));
    }
    Ok(name)
}

/// In-memory representation of an encoded sequence ready for disk persistence.
#[derive(Debug, Clone, PartialEq)]
pub struct EncodedSeq {
    /// Total sequence length in base pairs.
    pub sequence_length: u64,
    /// 2-bit packed canonical sequence bytes.
    pub packed_2bit: Vec<u8>,
    /// Ambiguity entries sorted by position: (position, IUPAC uppercase char).
    pub ambiguity_records: Vec<(u64, u8)>,
    /// True if Uracil was detected in the input stream.
    pub is_rna: bool,
}

/// On-disk store for 2-bit encoded chromosome sequences.
pub struct SeqStore {
    /// Root directory containing taxonomic/assembly `.seq` files.
    dir: PathBuf,
}

impl SeqStore {
    /// Opens (or creates) a [`SeqStore`] backed by `dir`.
    pub fn new(dir: &Path) -> Result<Self> {
        std::fs::create_dir_all(dir)?;
        Ok(Self {
            dir: dir.to_owned(),
        })
    }

    /// Returns the root directory path.
    pub fn root_dir(&self) -> &Path {
        &self.dir
    }

    /// Encodes `sequence` losslessly into `.seq v2` format and writes it atomically
    /// to `<dir>/<tax_id>/<assembly>/<chr>.seq`.
    pub fn store(&self, tax_id: u32, assembly: &str, chr: &str, sequence: &[u8]) -> Result<()> {
        let path = self.seq_path(tax_id, assembly, chr)?;
        let parent = path.parent().ok_or_else(|| {
            BioDbError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Invalid path parent",
            ))
        })?;
        std::fs::create_dir_all(parent)?;

        let encoded = Self::encode_v2(sequence)?;

        let packed_len = encoded.packed_2bit.len() as u64;
        let ambiguity_count = encoded.ambiguity_records.len() as u64;
        let ambiguity_len = ambiguity_count * 9;

        // Compute CRC32 over the payload bytes (packed + ambiguity)
        let mut crc_hasher_buf = Vec::with_capacity((packed_len + ambiguity_len) as usize);
        crc_hasher_buf.extend_from_slice(&encoded.packed_2bit);
        for &(pos, code) in &encoded.ambiguity_records {
            crc_hasher_buf.extend_from_slice(&pos.to_le_bytes());
            crc_hasher_buf.push(code);
        }
        let checksum = crc32(&crc_hasher_buf);

        // Build 64-byte v2 header
        let mut header = [0u8; HEADER_SIZE_V2];
        header[0..8].copy_from_slice(SEQ_V2_MAGIC);
        header[8..10].copy_from_slice(&SEQ_V2_VERSION.to_le_bytes());

        let mut flags = 0u16;
        if !encoded.ambiguity_records.is_empty() {
            flags |= 0x0001; // bit 0: contains ambiguity
        }
        if encoded.is_rna {
            flags |= 0x0002; // bit 1: RNA / Uracil present
        }
        header[10..12].copy_from_slice(&flags.to_le_bytes());
        header[12..20].copy_from_slice(&encoded.sequence_length.to_le_bytes());
        header[20..28].copy_from_slice(&packed_len.to_le_bytes());
        header[28..36].copy_from_slice(&ambiguity_count.to_le_bytes());
        header[36..44].copy_from_slice(&ambiguity_len.to_le_bytes());
        header[44..48].copy_from_slice(&checksum.to_le_bytes());
        // Remaining header[48..64] is zero-padded reserved

        // Atomic write via temporary file in same directory
        let temp_suffix = TEMP_FILE_COUNTER.fetch_add(1, Ordering::Relaxed);
        let temp_path = parent.join(format!(
            "{chr}.seq.tmp.{}_{temp_suffix}",
            std::process::id()
        ));

        {
            let temp_file = File::create(&temp_path)?;
            let mut writer = BufWriter::new(temp_file);
            writer.write_all(&header)?;
            writer.write_all(&encoded.packed_2bit)?;
            for &(pos, code) in &encoded.ambiguity_records {
                writer.write_all(&pos.to_le_bytes())?;
                writer.write_all(&[code])?;
            }
            writer.flush()?;
            let file = writer.into_inner().map_err(|e| e.into_error())?;
            file.sync_all()?;
        }

        // Atomic rename
        std::fs::rename(&temp_path, &path)?;
        Ok(())
    }

    /// Memory-maps the sequence file and returns decoded bases for the half-open interval `[start, end)`.
    pub fn fetch(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>> {
        let path = self.seq_path(tax_id, assembly, chr)?;
        if !path.exists() {
            return Err(BioDbError::NotFound(format!(
                "sequence file not found for tax_id {tax_id}, assembly {assembly}, chromosome {chr}"
            )));
        }
        let file = File::open(&path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        if mmap.len() < 8 {
            return Err(BioDbError::CorruptData(format!(
                "Sequence file for {chr} is too short to contain a valid header ({} bytes)",
                mmap.len()
            )));
        }

        if mmap.starts_with(SEQ_V2_MAGIC) {
            Self::fetch_v2(&mmap, chr, start, end)
        } else {
            Self::fetch_v1(&mmap, chr, start, end)
        }
    }

    /// Returns the total base length stored for chromosome `chr`.
    pub fn sequence_length(&self, tax_id: u32, assembly: &str, chr: &str) -> Result<u64> {
        let path = self.seq_path(tax_id, assembly, chr)?;
        if !path.exists() {
            return Err(BioDbError::NotFound(format!(
                "sequence file not found for tax_id {tax_id}, assembly {assembly}, chromosome {chr}"
            )));
        }
        let file = File::open(&path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        if mmap.len() < 8 {
            return Err(BioDbError::CorruptData(format!(
                "Sequence file for {chr} is too short to contain a valid header"
            )));
        }

        if mmap.starts_with(SEQ_V2_MAGIC) {
            if mmap.len() < HEADER_SIZE_V2 {
                return Err(BioDbError::CorruptData(format!(
                    "v2 header truncated for {chr}: expected at least {HEADER_SIZE_V2} bytes, got {}",
                    mmap.len()
                )));
            }
            let version = u16::from_le_bytes(mmap[8..10].try_into().unwrap());
            if version != SEQ_V2_VERSION {
                return Err(BioDbError::UnsupportedVersion(version));
            }
            Ok(u64::from_le_bytes(mmap[12..20].try_into().unwrap()))
        } else {
            Ok(u64::from_le_bytes(mmap[..8].try_into().unwrap()))
        }
    }

    /// Returns the `.seq` file path for chromosome `chr` within a given taxonomy & assembly,
    /// validating that path traversal is prevented.
    pub fn seq_path(&self, tax_id: u32, assembly: &str, chr: &str) -> Result<PathBuf> {
        let safe_assembly = sanitize_component(assembly, "assembly")?;
        let safe_chr = sanitize_component(chr, "chromosome")?;
        Ok(self
            .dir
            .join(format!("{tax_id}/{safe_assembly}/{safe_chr}.seq")))
    }

    /// Decodes a v2 sequence slice with full validation and ambiguity recovery.
    fn fetch_v2(mmap: &[u8], chr: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        if mmap.len() < HEADER_SIZE_V2 {
            return Err(BioDbError::CorruptData(format!(
                "v2 header truncated for {chr}: expected at least {HEADER_SIZE_V2} bytes, got {}",
                mmap.len()
            )));
        }

        let version = u16::from_le_bytes(mmap[8..10].try_into().unwrap());
        if version != SEQ_V2_VERSION {
            return Err(BioDbError::UnsupportedVersion(version));
        }

        let seq_len = u64::from_le_bytes(mmap[12..20].try_into().unwrap());
        let packed_len = u64::from_le_bytes(mmap[20..28].try_into().unwrap());
        let ambiguity_count = u64::from_le_bytes(mmap[28..36].try_into().unwrap());
        let ambiguity_len = u64::from_le_bytes(mmap[36..44].try_into().unwrap());
        let expected_checksum = u32::from_le_bytes(mmap[44..48].try_into().unwrap());

        // Validate payload boundaries
        let expected_packed = seq_len.div_ceil(4);
        if packed_len != expected_packed {
            return Err(BioDbError::CorruptData(format!(
                "Inconsistent packed length: header specifies {packed_len}, expected {expected_packed}"
            )));
        }

        let expected_ambiguity_len = ambiguity_count.saturating_mul(9);
        if ambiguity_len != expected_ambiguity_len {
            return Err(BioDbError::CorruptData(format!(
                "Inconsistent ambiguity length: header specifies {ambiguity_len}, expected {expected_ambiguity_len}"
            )));
        }

        let total_required = (HEADER_SIZE_V2 as u64)
            .checked_add(packed_len)
            .and_then(|v| v.checked_add(ambiguity_len))
            .ok_or_else(|| BioDbError::CorruptData("Payload size overflow".to_string()))?;

        if (mmap.len() as u64) < total_required {
            return Err(BioDbError::CorruptData(format!(
                "Truncated .seq file: expected {total_required} bytes, got {}",
                mmap.len()
            )));
        }

        // Checksum verification
        let payload_slice = &mmap[HEADER_SIZE_V2..total_required as usize];
        let computed_checksum = crc32(payload_slice);
        if computed_checksum != expected_checksum {
            return Err(BioDbError::ChecksumMismatch {
                expected: expected_checksum,
                found: computed_checksum,
            });
        }

        // Coordinate check
        if end > seq_len || start > end {
            return Err(BioDbError::InvalidCoordinate {
                chr: chr.to_owned(),
                start,
                end,
            });
        }

        // Slice 2-bit canonical payload
        let packed_slice = &mmap[HEADER_SIZE_V2..HEADER_SIZE_V2 + packed_len as usize];
        let mut out = Self::decode_2bit(packed_slice, start, end);

        // Apply ambiguity overrides within [start, end)
        if ambiguity_count > 0 {
            let ambiguity_start = HEADER_SIZE_V2 + packed_len as usize;
            let ambiguity_slice = &mmap[ambiguity_start..ambiguity_start + ambiguity_len as usize];

            let count = ambiguity_count as usize;
            // Binary search to find the lower bound index where position >= start
            let mut low = 0;
            let mut high = count;
            while low < high {
                let mid = (low + high) / 2;
                let pos =
                    u64::from_le_bytes(ambiguity_slice[mid * 9..mid * 9 + 8].try_into().unwrap());
                if pos < start {
                    low = mid + 1;
                } else {
                    high = mid;
                }
            }

            // Iterate until pos >= end
            for idx in low..count {
                let offset = idx * 9;
                let pos =
                    u64::from_le_bytes(ambiguity_slice[offset..offset + 8].try_into().unwrap());
                if pos >= end {
                    break;
                }
                let code = ambiguity_slice[offset + 8];
                let out_idx = (pos - start) as usize;
                out[out_idx] = code;
            }
        }

        Ok(out)
    }

    /// Decodes a v1 legacy sequence slice for backward compatibility.
    fn fetch_v1(mmap: &[u8], chr: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let total_bases = u64::from_le_bytes(mmap[..8].try_into().unwrap());
        let expected_payload = total_bases.div_ceil(4);
        let total_required = 8u64.saturating_add(expected_payload);

        if (mmap.len() as u64) < total_required {
            return Err(BioDbError::CorruptData(format!(
                "Truncated v1 .seq file: expected {total_required} bytes, got {}",
                mmap.len()
            )));
        }

        if end > total_bases || start > end {
            return Err(BioDbError::InvalidCoordinate {
                chr: chr.to_owned(),
                start,
                end,
            });
        }
        let data = &mmap[8..8 + expected_payload as usize];
        Ok(Self::decode_2bit(data, start, end))
    }

    /// Losslessly encodes an input sequence into 2-bit canonical bases and ambiguity sidecar.
    pub fn encode_v2(seq: &[u8]) -> Result<EncodedSeq> {
        let n = seq.len();
        let out_len = n.div_ceil(4);
        let mut packed_2bit = vec![0u8; out_len];
        let mut ambiguity_records = Vec::new();
        let mut is_rna = false;

        for (i, &raw_base) in seq.iter().enumerate() {
            let upper = raw_base.to_ascii_uppercase();
            let (bits, is_ambig) = match upper {
                b'A' => (0b00, false),
                b'C' => (0b01, false),
                b'G' => (0b10, false),
                b'T' => (0b11, false),
                b'U' => {
                    is_rna = true;
                    (0b11, true) // Canonical: Thymine bit pattern, recorded as U in ambiguity table
                }
                b'N' => (0b00, true),
                b'R' => (0b00, true), // A or G
                b'Y' => (0b01, true), // C or T
                b'S' => (0b10, true), // G or C
                b'W' => (0b11, true), // A or T
                b'K' => (0b10, true), // G or T
                b'M' => (0b00, true), // A or C
                b'B' => (0b01, true), // C, G or T
                b'D' => (0b00, true), // A, G or T
                b'H' => (0b00, true), // A, C or T
                b'V' => (0b00, true), // A, C or G
                invalid => {
                    return Err(BioDbError::Encoding(format!(
                        "Invalid IUPAC nucleotide '{c}' (0x{invalid:02X}) at position {i}",
                        c = invalid as char
                    )));
                }
            };

            let byte_idx = i / 4;
            let shift = 6 - (i % 4) * 2;
            packed_2bit[byte_idx] |= bits << shift;

            if is_ambig {
                ambiguity_records.push((i as u64, upper));
            }
        }

        Ok(EncodedSeq {
            sequence_length: n as u64,
            packed_2bit,
            ambiguity_records,
            is_rna,
        })
    }

    /// 2-bit-encodes a slice of DNA bytes for backward compatibility.
    pub fn encode(seq: &[u8]) -> Vec<u8> {
        Self::encode_v2(seq)
            .map(|enc| enc.packed_2bit)
            .unwrap_or_else(|_| {
                let n = seq.len();
                let out_len = n.div_ceil(4);
                let mut out = vec![0u8; out_len];
                for (i, &base) in seq.iter().enumerate() {
                    let bits: u8 = match base.to_ascii_uppercase() {
                        b'A' => 0b00,
                        b'C' => 0b01,
                        b'G' => 0b10,
                        b'T' => 0b11,
                        _ => 0b00,
                    };
                    let byte_idx = i / 4;
                    let shift = 6 - (i % 4) * 2;
                    out[byte_idx] |= bits << shift;
                }
                out
            })
    }

    /// Decodes 2-bit canonical data back to ASCII bytes for the half-open interval `[start, end)`.
    pub fn decode_2bit(data: &[u8], start: u64, end: u64) -> Vec<u8> {
        const BASES: [u8; 4] = *b"ACGT";
        let mut out = Vec::with_capacity((end - start) as usize);
        for pos in start..end {
            let byte_idx = (pos / 4) as usize;
            let shift = 6 - (pos % 4) * 2;
            let bits = (data[byte_idx] >> shift) & 0b11;
            out.push(BASES[bits as usize]);
        }
        out
    }

    /// Decodes 2-bit data back to ASCII bytes (legacy method signature).
    pub fn decode(data: &[u8], start: u64, end: u64) -> Vec<u8> {
        Self::decode_2bit(data, start, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_store_and_fetch_lossless_iupac() {
        let dir = tempdir().unwrap();
        let store = SeqStore::new(dir.path()).unwrap();

        // Sequence containing all IUPAC symbols and Uracil
        let seq = b"ACGTUNRYSWKMBDHVacgtunryswkmbdhv";
        store.store(9606, "GRCh38", "chr1", seq).unwrap();

        let fetched = store
            .fetch(9606, "GRCh38", "chr1", 0, seq.len() as u64)
            .unwrap();
        // Lowercase input is canonicalized to uppercase IUPAC
        let expected = b"ACGTUNRYSWKMBDHVACGTUNRYSWKMBDHV";
        assert_eq!(&fetched, expected);

        // Sub-slice verification
        let slice_u = store.fetch(9606, "GRCh38", "chr1", 4, 6).unwrap();
        assert_eq!(&slice_u, b"UN");

        let len = store.sequence_length(9606, "GRCh38", "chr1").unwrap();
        assert_eq!(len, seq.len() as u64);
    }

    #[test]
    fn test_invalid_base_rejected() {
        let bad_seq = b"ACGT123";
        let res = SeqStore::encode_v2(bad_seq);
        assert!(matches!(res, Err(BioDbError::Encoding(_))));
    }

    #[test]
    fn test_truncated_file_detection() {
        let dir = tempdir().unwrap();
        let store = SeqStore::new(dir.path()).unwrap();

        let seq = b"ACGTACGTACGT";
        store.store(9606, "GRCh38", "chr_trunc", seq).unwrap();

        let path = store.seq_path(9606, "GRCh38", "chr_trunc").unwrap();
        let data = std::fs::read(&path).unwrap();

        // Truncate file by removing trailing bytes
        let truncated = &data[..data.len() - 2];
        std::fs::write(&path, truncated).unwrap();

        let res = store.fetch(9606, "GRCh38", "chr_trunc", 0, 12);
        assert!(matches!(res, Err(BioDbError::CorruptData(_))));
    }

    #[test]
    fn test_path_traversal_protection() {
        let dir = tempdir().unwrap();
        let store = SeqStore::new(dir.path()).unwrap();

        let res = store.store(9606, "../evil", "chr1", b"ACGT");
        assert!(matches!(res, Err(BioDbError::PathTraversal(_))));

        let res2 = store.store(9606, "GRCh38", "sub/chr1", b"ACGT");
        assert!(matches!(res2, Err(BioDbError::PathTraversal(_))));
    }

    #[test]
    fn test_crc32_checksum_mismatch() {
        let dir = tempdir().unwrap();
        let store = SeqStore::new(dir.path()).unwrap();

        let seq = b"ACGTACGTACGT";
        store.store(9606, "GRCh38", "chr_crc", seq).unwrap();

        let path = store.seq_path(9606, "GRCh38", "chr_crc").unwrap();
        let mut data = std::fs::read(&path).unwrap();

        // Corrupt a byte in payload (after 64-byte header)
        let last_idx = data.len() - 1;
        data[last_idx] ^= 0xFF;
        std::fs::write(&path, data).unwrap();

        let res = store.fetch(9606, "GRCh38", "chr_crc", 0, 12);
        assert!(matches!(res, Err(BioDbError::ChecksumMismatch { .. })));
    }

    #[test]
    fn test_unsupported_version_detection() {
        let dir = tempdir().unwrap();
        let store = SeqStore::new(dir.path()).unwrap();

        let seq = b"ACGT";
        store.store(9606, "GRCh38", "chr_version", seq).unwrap();

        let path = store.seq_path(9606, "GRCh38", "chr_version").unwrap();
        let mut data = std::fs::read(&path).unwrap();

        // Change version from 2 to 99
        data[8] = 99;
        data[9] = 0;
        std::fs::write(&path, data).unwrap();

        let res = store.fetch(9606, "GRCh38", "chr_version", 0, 4);
        assert!(matches!(res, Err(BioDbError::UnsupportedVersion(99))));
    }

    #[test]
    fn test_bad_magic_detection() {
        let dir = tempdir().unwrap();
        let store = SeqStore::new(dir.path()).unwrap();

        let seq = b"ACGT";
        store.store(9606, "GRCh38", "chr_magic", seq).unwrap();

        let path = store.seq_path(9606, "GRCh38", "chr_magic").unwrap();
        let mut data = std::fs::read(&path).unwrap();

        // Corrupt magic bytes
        data[0] = b'X';
        data[1] = b'Y';
        std::fs::write(&path, data).unwrap();

        let res = store.fetch(9606, "GRCh38", "chr_magic", 0, 4);
        assert!(matches!(res, Err(BioDbError::CorruptData(_))));
    }
}
