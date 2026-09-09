//! K-mer hash index for fast exact-match sequence lookup.
//!
//! Each k-mer is mapped to a list of `(chromosome, position)` pairs so that
//! a query k-mer can be located in O(1) expected time.

use std::collections::HashMap;

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
