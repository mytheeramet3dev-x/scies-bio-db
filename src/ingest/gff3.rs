//! GFF3 annotation ingestion pipeline.
//!
//! Parses GFF3 records and inserts gene, transcript, and exon rows into the
//! [`SqliteStore`].

use std::path::Path;

use crate::storage::SqliteStore;
use crate::Result;

/// Ingests genomic feature annotations from GFF3 files into the [`SqliteStore`].
pub struct Gff3Ingestor;

impl Gff3Ingestor {
    /// Creates a new [`Gff3Ingestor`].
    pub fn new() -> Self {
        Self
    }

    /// Parses `path` as a GFF3 file and inserts all recognised features into
    /// `store`.
    ///
    /// Returns the total number of features successfully inserted.
    pub fn ingest_file(&self, path: &Path, store: &SqliteStore) -> Result<u64> {
        let _ = (path, store);
        todo!("parse GFF3 lines, build gene/transcript/exon structs, insert via store.conn")
    }
}

impl Default for Gff3Ingestor {
    fn default() -> Self {
        Self::new()
    }
}
