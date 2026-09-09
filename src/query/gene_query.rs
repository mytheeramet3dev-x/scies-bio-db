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

    /// Finds a gene by its HGNC symbol across all organisms.
    pub fn find_by_symbol(&self, symbol: &str) -> Result<Gene> {
        self.sql_store
            .find_gene_by_symbol_any(symbol)?
            .ok_or_else(|| BioDbError::NotFound(format!("gene with symbol '{symbol}' not found")))
    }

    /// Finds a gene by tax_id and HGNC symbol (species-isolated).
    pub fn find_by_symbol_scoped(&self, tax_id: u32, symbol: &str) -> Result<Gene> {
        self.sql_store
            .find_gene_by_symbol(tax_id, symbol)?
            .ok_or_else(|| {
                BioDbError::NotFound(format!(
                    "gene with symbol '{symbol}' not found for tax_id {tax_id}"
                ))
            })
    }

    /// Finds a gene by its stable Ensembl or HGNC gene identifier across all assemblies.
    pub fn find_by_id(&self, gene_id: &str) -> Result<Gene> {
        self.sql_store
            .find_gene_by_id(gene_id)?
            .ok_or_else(|| BioDbError::NotFound(format!("gene with ID '{gene_id}' not found")))
    }

    /// Finds a gene by its stable identifier scoped to a specific assembly.
    pub fn find_by_id_scoped(&self, gene_id: &str, assembly: &str) -> Result<Gene> {
        self.sql_store
            .find_gene_by_id_scoped(gene_id, assembly)?
            .ok_or_else(|| {
                BioDbError::NotFound(format!(
                    "gene with ID '{gene_id}' not found in assembly '{assembly}'"
                ))
            })
    }

    /// Returns all genes whose locus overlaps interval `[start, end)` on chromosome `chr`,
    /// strictly scoped to `tax_id` and `assembly`.
    pub fn find_by_region(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<Gene>> {
        self.sql_store
            .find_genes_by_region(tax_id, assembly, chr, start, end)
    }

    /// Unscoped region query (convenience fallback).
    pub fn find_by_region_unscoped(&self, chr: &str, start: u64, end: u64) -> Result<Vec<Gene>> {
        self.sql_store
            .find_genes_by_region_unscoped(chr, start, end)
    }

    /// Returns all transcripts belonging to gene `gene_id`.
    pub fn transcripts_of(&self, gene_id: &str) -> Result<Vec<Transcript>> {
        self.sql_store.find_transcripts_by_gene(gene_id)
    }

    /// Returns all transcripts belonging to gene `gene_id` within a specific `assembly`.
    pub fn transcripts_of_scoped(&self, gene_id: &str, assembly: &str) -> Result<Vec<Transcript>> {
        self.sql_store
            .find_transcripts_by_gene_scoped(gene_id, assembly)
    }
}
