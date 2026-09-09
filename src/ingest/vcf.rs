//! VCF variant ingestion pipeline.
//!
//! Parses VCF records and inserts [`VariantRecord`] rows into the
//! [`SqliteStore`].

use std::path::Path;

use crate::storage::SqliteStore;
use crate::Result;

/// Ingests genomic variants from VCF files into the [`SqliteStore`].
pub struct VcfIngestor;

impl VcfIngestor {
    /// Creates a new [`VcfIngestor`].
    pub fn new() -> Self {
        Self
    }

    /// Parses `path` as a VCF file and inserts all variant records into `store`.
    ///
    /// Returns the total number of variants successfully inserted.
    pub fn ingest_file(&self, path: &Path, store: &SqliteStore) -> Result<u64> {
        let _ = (path, store);
        todo!("parse VCF header + data lines, build VariantRecord, insert via store.conn")
    }
}

impl Default for VcfIngestor {
    fn default() -> Self {
        Self::new()
    }
}
