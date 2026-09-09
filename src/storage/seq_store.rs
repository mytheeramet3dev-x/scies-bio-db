//! 2-bit encoded genomic sequence storage organized by taxonomy and assembly.
//!
//! Each base is packed into 2 bits:
//!
//! | Base | Bits |
//! |------|------|
//! | A    | 00   |
//! | C    | 01   |
//! | G    | 10   |
//! | T    | 11   |
//!
//! Sequences are stored as flat `.seq` files under:
//! `<dir>/<tax_id>/<assembly>/<chr>.seq`

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use memmap2::MmapOptions;

use crate::{BioDbError, Result};

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

    /// Encodes `sequence` in 2-bit format and writes it to `<dir>/<tax_id>/<assembly>/<chr>.seq`.
    pub fn store(&self, tax_id: u32, assembly: &str, chr: &str, sequence: &[u8]) -> Result<()> {
        let path = self.seq_path(tax_id, assembly, chr);
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let file = std::fs::File::create(&path)?;
        let mut writer = BufWriter::new(file);
        let len = sequence.len() as u64;
        writer.write_all(&len.to_le_bytes())?;
        let encoded = Self::encode(sequence);
        writer.write_all(&encoded)?;
        writer.flush()?;
        Ok(())
    }

    /// Memory-maps `<dir>/<tax_id>/<assembly>/<chr>.seq` and returns the decoded bases for the
    /// half-open interval `[start, end)`.
    pub fn fetch(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>> {
        let path = self.seq_path(tax_id, assembly, chr);
        if !path.exists() {
            return Err(BioDbError::NotFound(format!(
                "sequence file not found for tax_id {tax_id}, assembly {assembly}, chromosome {chr}"
            )));
        }
        let file = std::fs::File::open(&path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        if mmap.len() < 8 {
            return Err(BioDbError::Encoding(format!(
                "seq file for {chr} is too short to contain a header"
            )));
        }
        let total_bases = u64::from_le_bytes(mmap[..8].try_into().unwrap());
        if end > total_bases || start > end {
            return Err(BioDbError::InvalidCoordinate {
                chr: chr.to_owned(),
                start,
                end,
            });
        }
        let data = &mmap[8..];
        Ok(Self::decode(data, start, end))
    }

    /// Returns the total base length stored for chromosome `chr`.
    pub fn sequence_length(&self, tax_id: u32, assembly: &str, chr: &str) -> Result<u64> {
        let path = self.seq_path(tax_id, assembly, chr);
        if !path.exists() {
            return Err(BioDbError::NotFound(format!(
                "sequence file not found for tax_id {tax_id}, assembly {assembly}, chromosome {chr}"
            )));
        }
        let file = std::fs::File::open(&path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        if mmap.len() < 8 {
            return Err(BioDbError::Encoding(format!(
                "seq file for {chr} is too short to contain a header"
            )));
        }
        Ok(u64::from_le_bytes(mmap[..8].try_into().unwrap()))
    }

    /// Returns the `.seq` file path for chromosome `chr` within a given taxonomy & assembly.
    pub fn seq_path(&self, tax_id: u32, assembly: &str, chr: &str) -> PathBuf {
        self.dir.join(format!("{tax_id}/{assembly}/{chr}.seq"))
    }

    /// 2-bit-encodes a slice of ASCII DNA bytes (`A`, `C`, `G`, `T`).
    pub fn encode(seq: &[u8]) -> Vec<u8> {
        let n = seq.len();
        let out_len = n.div_ceil(4);
        let mut out = vec![0u8; out_len];
        for (i, &base) in seq.iter().enumerate() {
            let bits: u8 = match base.to_ascii_uppercase() {
                b'A' => 0b00,
                b'C' => 0b01,
                b'G' => 0b10,
                b'T' => 0b11,
                _ => 0b00, // treat ambiguous bases as A
            };
            let byte_idx = i / 4;
            let shift = 6 - (i % 4) * 2; // MSB-first packing
            out[byte_idx] |= bits << shift;
        }
        out
    }

    /// Decodes 2-bit data back to ASCII bytes for the half-open interval `[start, end)`.
    pub fn decode(data: &[u8], start: u64, end: u64) -> Vec<u8> {
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_store_and_fetch() {
        let dir = tempdir().unwrap();
        let store = SeqStore::new(dir.path()).unwrap();

        let seq = b"ACGTACGTNNACGT";
        store.store(9606, "GRCh38", "chr1", seq).unwrap();

        let fetched = store.fetch(9606, "GRCh38", "chr1", 0, 8).unwrap();
        assert_eq!(&fetched, b"ACGTACGT");

        let fetched_slice = store.fetch(9606, "GRCh38", "chr1", 4, 8).unwrap();
        assert_eq!(&fetched_slice, b"ACGT");
    }
}
