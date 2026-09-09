//! Variant and clinical association queries.

use std::sync::Arc;

use crate::model::variant::{ClinicalVariant, VariantRecord};
use crate::storage::SqliteStore;
use crate::Result;

/// Handles queries against the variant records and clinical variants tables.
pub struct VariantQuery {
    #[allow(dead_code)]
    pub(crate) sql_store: Arc<SqliteStore>,
}

impl VariantQuery {
    /// Creates a new [`VariantQuery`] backed by `sql_store`.
    pub fn new(sql_store: Arc<SqliteStore>) -> Self {
        Self { sql_store }
    }

    /// Finds a [`VariantRecord`] by its dbSNP rs identifier.
    pub fn find_by_id(&self, variant_id: &str) -> Result<VariantRecord> {
        let _ = variant_id;
        todo!("SELECT * FROM variant_records WHERE variant_id = ?1")
    }

    /// Returns all variants located within the half-open interval
    /// `[start, end)` on chromosome `chr`.
    pub fn find_by_region(&self, chr: &str, start: u64, end: u64) -> Result<Vec<VariantRecord>> {
        let _ = (chr, start, end);
        todo!("SELECT * FROM variant_records WHERE chr = ?1 AND position >= ?2 AND position < ?3")
    }

    /// Returns all clinical variant interpretations for a given gene symbol.
    pub fn find_by_gene(&self, gene_symbol: &str) -> Result<Vec<ClinicalVariant>> {
        let _ = gene_symbol;
        todo!("SELECT * FROM clinical_variants WHERE gene_symbol = ?1")
    }
}
