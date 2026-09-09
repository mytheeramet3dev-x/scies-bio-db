//! Protein entry and domain queries.

use std::sync::Arc;

use crate::model::protein::{Domain, ProteinEntry};
use crate::storage::SqliteStore;
use crate::Result;

/// Handles queries against the protein entries and domains tables.
pub struct ProteinQuery {
    #[allow(dead_code)]
    pub(crate) sql_store: Arc<SqliteStore>,
}

impl ProteinQuery {
    /// Creates a new [`ProteinQuery`] backed by `sql_store`.
    pub fn new(sql_store: Arc<SqliteStore>) -> Self {
        Self { sql_store }
    }

    /// Finds a [`ProteinEntry`] by its primary UniProtKB accession.
    pub fn find_by_accession(&self, accession: &str) -> Result<ProteinEntry> {
        let _ = accession;
        todo!("SELECT * FROM protein_entries WHERE uniprot_accession = ?1")
    }

    /// Returns all [`ProteinEntry`] records linked to a given gene symbol.
    pub fn find_by_gene(&self, gene_symbol: &str) -> Result<Vec<ProteinEntry>> {
        let _ = gene_symbol;
        todo!("SELECT * FROM protein_entries WHERE gene_symbol = ?1")
    }

    /// Returns all [`Domain`] records annotated on the protein with the given
    /// UniProtKB accession.
    pub fn domains_of(&self, accession: &str) -> Result<Vec<Domain>> {
        let _ = accession;
        todo!("SELECT * FROM protein_domains WHERE uniprot_accession = ?1")
    }
}
