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
        self.sql_store
            .find_variant_by_id(variant_id)?
            .ok_or_else(|| {
                crate::BioDbError::NotFound(format!("variant with ID '{variant_id}' not found"))
            })
    }

    /// Finds a [`VariantRecord`] by its dbSNP rs identifier scoped to an assembly.
    pub fn find_by_id_scoped(&self, variant_id: &str, assembly: &str) -> Result<VariantRecord> {
        self.sql_store
            .find_variant_by_id_scoped(variant_id, assembly)?
            .ok_or_else(|| {
                crate::BioDbError::NotFound(format!(
                    "variant with ID '{variant_id}' not found in assembly '{assembly}'"
                ))
            })
    }

    /// Returns all variants located within interval `[start, end)` on chromosome `chr`,
    /// scoped to `tax_id` and `assembly`.
    pub fn find_by_region(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<VariantRecord>> {
        self.sql_store
            .find_variants_by_region(tax_id, assembly, chr, start, end)
    }

    /// Returns all clinical variant interpretations for a given gene symbol.
    pub fn find_by_gene(&self, gene_symbol: &str) -> Result<Vec<ClinicalVariant>> {
        self.sql_store.find_clinical_variants_by_gene(gene_symbol)
    }
}
