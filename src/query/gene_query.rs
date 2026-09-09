//! Gene and transcript queries.

use std::sync::Arc;

use crate::model::gene::{Gene, Transcript};
use crate::storage::SqliteStore;
use crate::{BioDbError, Result};

/// Handles queries against the genes and transcripts tables.
pub struct GeneQuery {
    pub(crate) sql_store: Arc<SqliteStore>,
}

impl GeneQuery {
    /// Creates a new [`GeneQuery`] backed by `sql_store`.
    pub fn new(sql_store: Arc<SqliteStore>) -> Self {
        Self { sql_store }
    }

    /// Finds a gene by its HGNC symbol (case-sensitive).
    pub fn find_by_symbol(&self, symbol: &str) -> Result<Gene> {
        self.sql_store
            .find_gene_by_symbol(symbol)?
            .ok_or_else(|| BioDbError::NotFound(format!("gene with symbol '{symbol}' not found")))
    }

    /// Finds a gene by its stable Ensembl or HGNC gene identifier.
    pub fn find_by_id(&self, gene_id: &str) -> Result<Gene> {
        self.sql_store
            .find_gene_by_id(gene_id)?
            .ok_or_else(|| BioDbError::NotFound(format!("gene with ID '{gene_id}' not found")))
    }

    /// Returns all genes whose locus overlaps the half-open interval
    /// `[start, end)` on chromosome `chr`.
    pub fn find_by_region(&self, chr: &str, start: u64, end: u64) -> Result<Vec<Gene>> {
        self.sql_store.find_genes_by_region(chr, start, end)
    }

    /// Returns all transcripts belonging to gene `gene_id`.
    pub fn transcripts_of(&self, gene_id: &str) -> Result<Vec<Transcript>> {
        self.sql_store.find_transcripts_by_gene(gene_id)
    }
}
