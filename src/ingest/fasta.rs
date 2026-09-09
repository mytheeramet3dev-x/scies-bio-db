//! FASTA ingestion pipeline.
//!
//! Reads one or more FASTA records, 2-bit-encodes each sequence, and stores
//! the result in a [`SeqStore`].

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use scies_bio_th::io::FastaReader;

use crate::model::genome::Chromosome;
use crate::storage::{SeqStore, SqliteStore};
use crate::Result;

/// Ingests chromosome sequences from FASTA files into the [`SeqStore`].
pub struct FastaIngestor {
    /// Shared reference to the sequence storage backend.
    seq_store: Arc<SeqStore>,
}

impl FastaIngestor {
    /// Creates a new [`FastaIngestor`] that will write sequences to `seq_store`.
    pub fn new(seq_store: Arc<SeqStore>) -> Self {
        Self { seq_store }
    }

    /// Parses `path` as a FASTA file and stores every sequence found under `(tax_id, assembly)`.
    ///
    /// If `sql_store` is provided, also inserts records into the `chromosomes` table.
    /// Returns the list of chromosome names that were successfully ingested.
    pub fn ingest_file(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
        sql_store: Option<&SqliteStore>,
    ) -> Result<Vec<String>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let fasta_reader = FastaReader::new(reader)?;

        let mut ingested_chromosomes = Vec::new();

        for item in fasta_reader {
            let record = item?;
            // Use first word of the header as chromosome name (e.g. "chr1" from ">chr1 Homo sapiens...")
            let chr_name = record
                .id
                .split_whitespace()
                .next()
                .unwrap_or(&record.id)
                .to_string();

            let seq_len = record.sequence.len() as u64;
            self.seq_store
                .store(tax_id, assembly, &chr_name, &record.sequence)?;

            if let Some(store) = sql_store {
                let is_mito = chr_name.eq_ignore_ascii_case("chrM")
                    || chr_name.eq_ignore_ascii_case("MT")
                    || chr_name.eq_ignore_ascii_case("chrMT");

                let chrom = Chromosome {
                    name: chr_name.clone(),
                    length: seq_len,
                    assembly_name: assembly.to_string(),
                    is_mitochondrial: is_mito,
                };
                store.insert_chromosome(&chrom)?;
            }

            ingested_chromosomes.push(chr_name);
        }

        Ok(ingested_chromosomes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_fasta_ingestion() {
        let dir = tempdir().unwrap();
        let seq_store = Arc::new(SeqStore::new(dir.path()).unwrap());
        let ingestor = FastaIngestor::new(seq_store.clone());

        let mut fasta_file = NamedTempFile::new().unwrap();
        writeln!(fasta_file, ">chr1 Human Chromosome 1").unwrap();
        writeln!(fasta_file, "ACGTACGT").unwrap();
        writeln!(fasta_file, ">chr2 Human Chromosome 2").unwrap();
        writeln!(fasta_file, "TGCATGCA").unwrap();
        fasta_file.flush().unwrap();

        let ingested = ingestor
            .ingest_file(fasta_file.path(), 9606, "GRCh38", None)
            .unwrap();

        assert_eq!(ingested, vec!["chr1", "chr2"]);

        let fetched1 = seq_store.fetch(9606, "GRCh38", "chr1", 0, 8).unwrap();
        assert_eq!(&fetched1, b"ACGTACGT");

        let fetched2 = seq_store.fetch(9606, "GRCh38", "chr2", 0, 8).unwrap();
        assert_eq!(&fetched2, b"TGCATGCA");
    }
}
