//! Genome-level queries: sequence fetch and chromosome listing.

use std::sync::Arc;

use crate::index::KmerIndex;
use crate::model::genome::Chromosome;
use crate::storage::{SeqStore, SqliteStore};
use crate::{BioDbError, Result};

/// Handles queries against chromosome sequences and assembly metadata.
pub struct GenomeQuery {
    /// Shared sequence storage backend.
    pub(crate) seq_store: Arc<SeqStore>,
    /// Shared SQLite metadata store.
    pub(crate) sql_store: Arc<SqliteStore>,
}

impl GenomeQuery {
    /// Creates a new [`GenomeQuery`].
    pub fn new(seq_store: Arc<SeqStore>, sql_store: Arc<SqliteStore>) -> Self {
        Self {
            seq_store,
            sql_store,
        }
    }

    /// Fetches the raw DNA sequence for `chr` over the half-open interval
    /// `[start, end)` and returns it as a [`scies_bio_th::Dna`] value.
    pub fn fetch_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<scies_bio_th::Dna> {
        let raw_bytes = self.seq_store.fetch(tax_id, assembly, chr, start, end)?;
        scies_bio_th::Dna::new(&raw_bytes).map_err(BioDbError::Bio)
    }

    /// Fetches raw sequence bytes for `chr` over the half-open interval `[start, end)`.
    ///
    /// Preserves full IUPAC degenerate symbols and Uracil without converting to [`scies_bio_th::Dna`].
    pub fn fetch_raw_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>> {
        self.seq_store.fetch(tax_id, assembly, chr, start, end)
    }

    /// Fetches the full DNA sequence for `chr` and returns it as a [`scies_bio_th::Dna`] value.
    pub fn fetch_full_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
    ) -> Result<scies_bio_th::Dna> {
        let len = self.seq_store.sequence_length(tax_id, assembly, chr)?;
        self.fetch_sequence(tax_id, assembly, chr, 0, len)
    }

    /// Fetches the complete raw sequence bytes for `chr`.
    pub fn fetch_full_raw_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
    ) -> Result<Vec<u8>> {
        let len = self.seq_store.sequence_length(tax_id, assembly, chr)?;
        self.fetch_raw_sequence(tax_id, assembly, chr, 0, len)
    }

    /// Returns metadata for every chromosome in assembly `assembly`.
    pub fn list_chromosomes(&self, assembly: &str) -> Result<Vec<Chromosome>> {
        self.sql_store.list_chromosomes(assembly)
    }

    /// Builds an in-memory [`KmerIndex`] of length `k` for chromosome `chr`.
    pub fn build_kmer_index(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        k: usize,
    ) -> Result<KmerIndex> {
        let raw = self.fetch_full_raw_sequence(tax_id, assembly, chr)?;
        let mut index = KmerIndex::new(k);
        index.index_sequence(chr, &raw);
        Ok(index)
    }

    /// Computes the Jaccard k-mer similarity between chromosome `chr` and `query_sequence`.
    pub fn estimate_similarity(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        query_sequence: &[u8],
        k: usize,
    ) -> Result<f64> {
        let index = self.build_kmer_index(tax_id, assembly, chr, k)?;
        Ok(index.jaccard_similarity(query_sequence))
    }
}
