//! K-mer hash index for fast exact-match sequence lookup.
//!
//! Each k-mer is mapped to a list of `(chromosome, position)` pairs so that
//! a query k-mer can be located in O(1) expected time.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use crate::storage::seq_store::crc32;
use crate::{BioDbError, Result};

/// Magic header bytes for persistent `.kdx` k-mer index files.
pub const KDX_MAGIC: &[u8; 8] = b"SCIESKDX";
/// Format version for `.kdx` files.
pub const KDX_VERSION: u16 = 1;

/// A located match of a k-mer sequence hit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerHit {
    /// Chromosome name where hit occurs.
    pub chr: String,
    /// 0-based start coordinate.
    pub position: u64,
    /// Count of shared k-mers in the matching window.
    pub hit_count: usize,
}

/// K-mer hash index.
///
/// The table maps a 64-bit rolling hash of each k-mer to the list of genomic
/// positions where that k-mer occurs.
pub struct KmerIndex {
    /// K-mer length.
    pub k: usize,
    /// Hash → list of `(chr, 0-based start position)` hit locations.
    table: HashMap<u64, Vec<(String, u64)>>,
}

impl KmerIndex {
    /// Creates a new [`KmerIndex`] for k-mers of length `k`.
    pub fn new(k: usize) -> Self {
        Self {
            k,
            table: HashMap::new(),
        }
    }

    /// Indexes all k-mers in `sequence` from chromosome `chr`.
    ///
    /// Positions are 0-based and refer to the start of each k-mer window.
    pub fn index_sequence(&mut self, chr: &str, sequence: &[u8]) {
        if sequence.len() < self.k {
            return;
        }
        for i in 0..=(sequence.len() - self.k) {
            let kmer = &sequence[i..i + self.k];
            let h = Self::hash_kmer(kmer);
            self.table
                .entry(h)
                .or_default()
                .push((chr.to_owned(), i as u64));
        }
    }

    /// Looks up all genomic positions for `kmer`.
    ///
    /// Returns `None` when the k-mer is absent from the index.
    pub fn lookup(&self, kmer: &[u8]) -> Option<&Vec<(String, u64)>> {
        let h = Self::hash_kmer(kmer);
        self.table.get(&h)
    }

    /// Computes the Jaccard similarity coefficient between the indexed k-mers and `query_sequence`:
    /// $$J(A, B) = \frac{|A \cap B|}{|A \cup B|}$$
    ///
    /// Returns a float in `[0.0, 1.0]`. Returns `0.0` if either set is empty.
    pub fn jaccard_similarity(&self, query_sequence: &[u8]) -> f64 {
        if query_sequence.len() < self.k || self.table.is_empty() {
            return 0.0;
        }

        let mut query_set = HashSet::new();
        for i in 0..=(query_sequence.len() - self.k) {
            let kmer = &query_sequence[i..i + self.k];
            query_set.insert(Self::hash_kmer(kmer));
        }

        let mut intersection = 0usize;
        for h in &query_set {
            if self.table.contains_key(h) {
                intersection += 1;
            }
        }

        let union = self.table.len() + query_set.len() - intersection;
        if union == 0 {
            0.0
        } else {
            intersection as f64 / union as f64
        }
    }

    /// Counts how many k-mers in `query_sequence` match at least once in the index.
    pub fn count_shared_kmers(&self, query_sequence: &[u8]) -> usize {
        if query_sequence.len() < self.k {
            return 0;
        }
        let mut count = 0;
        for i in 0..=(query_sequence.len() - self.k) {
            let kmer = &query_sequence[i..i + self.k];
            let h = Self::hash_kmer(kmer);
            if self.table.contains_key(&h) {
                count += 1;
            }
        }
        count
    }

    /// Locates genomic regions with the highest density of matching k-mers.
    pub fn find_matches(&self, query_sequence: &[u8], min_hits: usize) -> Vec<KmerHit> {
        if query_sequence.len() < self.k {
            return Vec::new();
        }

        let mut hit_counts: HashMap<(String, u64), usize> = HashMap::new();
        for i in 0..=(query_sequence.len() - self.k) {
            let kmer = &query_sequence[i..i + self.k];
            let h = Self::hash_kmer(kmer);
            if let Some(positions) = self.table.get(&h) {
                for (chr, pos) in positions {
                    // Cluster by rough 100bp windows to identify candidate homologous blocks
                    let window = pos / 100 * 100;
                    *hit_counts.entry((chr.clone(), window)).or_insert(0) += 1;
                }
            }
        }

        let mut matches: Vec<KmerHit> = hit_counts
            .into_iter()
            .filter(|(_, count)| *count >= min_hits)
            .map(|((chr, position), hit_count)| KmerHit {
                chr,
                position,
                hit_count,
            })
            .collect();

        matches.sort_by_key(|a| std::cmp::Reverse(a.hit_count));
        matches
    }

    /// Computes a collision-free 2-bit packed `u64` key for canonical DNA k-mers ($k \le 31$),
    /// falling back to a polynomial rolling hash if $k > 31$ or non-canonical bases are present.
    fn hash_kmer(kmer: &[u8]) -> u64 {
        if kmer.len() <= 31 {
            let mut packed: u64 = 1u64 << 62; // Tag bit distinguishing 2-bit lossless keys
            let mut all_canonical = true;
            for &byte in kmer {
                let bits = match byte.to_ascii_uppercase() {
                    b'A' => 0b00u64,
                    b'C' => 0b01u64,
                    b'G' => 0b10u64,
                    b'T' | b'U' => 0b11u64,
                    _ => {
                        all_canonical = false;
                        break;
                    }
                };
                packed = (packed << 2) | bits;
            }
            if all_canonical {
                return packed;
            }
        }

        const BASE: u64 = 131;
        let mut h: u64 = 0;
        for &byte in kmer {
            h = h.wrapping_mul(BASE).wrapping_add(byte as u64);
        }
        h & !(1u64 << 62)
    }

    /// Persists the [`KmerIndex`] to a binary `.kdx` file on disk with CRC32 integrity verification.
    pub fn save_to_file(&self, path: &Path) -> Result<()> {
        let mut payload = Vec::new();
        payload.extend_from_slice(&(self.k as u32).to_le_bytes());
        payload.extend_from_slice(&(self.table.len() as u64).to_le_bytes());

        // Sort keys for deterministic binary output
        let mut keys: Vec<u64> = self.table.keys().copied().collect();
        keys.sort_unstable();

        for key in keys {
            let hits = &self.table[&key];
            payload.extend_from_slice(&key.to_le_bytes());
            payload.extend_from_slice(&(hits.len() as u32).to_le_bytes());
            for (chr, pos) in hits {
                let chr_bytes = chr.as_bytes();
                payload.extend_from_slice(&(chr_bytes.len() as u16).to_le_bytes());
                payload.extend_from_slice(chr_bytes);
                payload.extend_from_slice(&pos.to_le_bytes());
            }
        }

        let checksum = crc32(&payload);
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writer.write_all(KDX_MAGIC)?;
        writer.write_all(&KDX_VERSION.to_le_bytes())?;
        writer.write_all(&checksum.to_le_bytes())?;
        writer.write_all(&payload)?;
        writer.flush()?;
        Ok(())
    }

    /// Loads a [`KmerIndex`] from a binary `.kdx` file on disk, verifying magic bytes and CRC32.
    pub fn load_from_file(path: &Path) -> Result<Self> {
        let mut file = File::open(path)?;
        let mut buf = Vec::new();
        file.read_to_end(&mut buf)?;

        if buf.len() < 26 {
            return Err(BioDbError::CorruptData(
                "Truncated .kdx file header".to_string(),
            ));
        }
        if &buf[0..8] != KDX_MAGIC {
            return Err(BioDbError::CorruptData(
                "Invalid .kdx magic bytes".to_string(),
            ));
        }
        let version = u16::from_le_bytes(buf[8..10].try_into().unwrap());
        if version != KDX_VERSION {
            return Err(BioDbError::UnsupportedVersion(version));
        }
        let expected_crc = u32::from_le_bytes(buf[10..14].try_into().unwrap());
        let payload = &buf[14..];
        let actual_crc = crc32(payload);
        if actual_crc != expected_crc {
            return Err(BioDbError::ChecksumMismatch {
                expected: expected_crc,
                found: actual_crc,
            });
        }

        let k = u32::from_le_bytes(payload[0..4].try_into().unwrap()) as usize;
        let bucket_count = u64::from_le_bytes(payload[4..12].try_into().unwrap()) as usize;
        let mut offset = 12usize;
        let mut table = HashMap::with_capacity(bucket_count);

        for _ in 0..bucket_count {
            if offset + 12 > payload.len() {
                return Err(BioDbError::CorruptData(
                    "Unexpected EOF reading .kdx bucket".to_string(),
                ));
            }
            let key = u64::from_le_bytes(payload[offset..offset + 8].try_into().unwrap());
            let hit_count =
                u32::from_le_bytes(payload[offset + 8..offset + 12].try_into().unwrap()) as usize;
            offset += 12;

            let mut hits = Vec::with_capacity(hit_count);
            for _ in 0..hit_count {
                if offset + 2 > payload.len() {
                    return Err(BioDbError::CorruptData(
                        "Unexpected EOF reading .kdx chr length".to_string(),
                    ));
                }
                let chr_len =
                    u16::from_le_bytes(payload[offset..offset + 2].try_into().unwrap()) as usize;
                offset += 2;
                if offset + chr_len + 8 > payload.len() {
                    return Err(BioDbError::CorruptData(
                        "Unexpected EOF reading .kdx hit record".to_string(),
                    ));
                }
                let chr = String::from_utf8_lossy(&payload[offset..offset + chr_len]).to_string();
                offset += chr_len;
                let pos = u64::from_le_bytes(payload[offset..offset + 8].try_into().unwrap());
                offset += 8;
                hits.push((chr, pos));
            }
            table.insert(key, hits);
        }

        Ok(Self { k, table })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_index() {
        let mut index = KmerIndex::new(4);
        let seq = b"ACGTACGT";
        index.index_sequence("chr1", seq);

        let hits = index.lookup(b"ACGT").unwrap();
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0], ("chr1".to_string(), 0));
        assert_eq!(hits[1], ("chr1".to_string(), 4));

        let hits_cgt = index.lookup(b"CGTA").unwrap();
        assert_eq!(hits_cgt.len(), 1);
        assert_eq!(hits_cgt[0], ("chr1".to_string(), 1));

        assert!(index.lookup(b"TTTT").is_none());

        // Test similarity
        let identical_score = index.jaccard_similarity(b"ACGTACGT");
        assert!((identical_score - 1.0).abs() < 1e-6);

        let shared_count = index.count_shared_kmers(b"ACGTGGGG");
        assert_eq!(shared_count, 1); // only ACGT matches

        let matches = index.find_matches(b"ACGTACGT", 1);
        assert!(!matches.is_empty());
        assert_eq!(matches[0].chr, "chr1");

        // Test .kdx persistence roundtrip
        let tmp = tempfile::NamedTempFile::new().unwrap();
        index.save_to_file(tmp.path()).unwrap();
        let loaded = KmerIndex::load_from_file(tmp.path()).unwrap();
        assert_eq!(loaded.k, 4);
        assert_eq!(loaded.lookup(b"ACGT").unwrap().len(), 2);
    }
}
