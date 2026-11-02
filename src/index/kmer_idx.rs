//! K-mer hash index for fast exact-match sequence lookup.
//!
//! Each k-mer is mapped to a list of `(chromosome, position)` pairs so that
//! a query k-mer can be located in O(1) expected time.

use std::collections::{HashMap, HashSet};

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

    /// Computes a simple polynomial rolling hash for `kmer`.
    ///
    /// Uses base 131 and modulus 2^64 (wrapping arithmetic).
    fn hash_kmer(kmer: &[u8]) -> u64 {
        const BASE: u64 = 131;
        let mut h: u64 = 0;
        for &byte in kmer {
            h = h.wrapping_mul(BASE).wrapping_add(byte as u64);
        }
        h
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
    }
}
