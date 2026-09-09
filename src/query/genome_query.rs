//! Genome-level queries: sequence fetch and chromosome listing.

use std::sync::Arc;

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

    /// Returns metadata for every chromosome in assembly `assembly`.
    pub fn list_chromosomes(&self, assembly: &str) -> Result<Vec<Chromosome>> {
        self.sql_store.list_chromosomes(assembly)
    }
}
