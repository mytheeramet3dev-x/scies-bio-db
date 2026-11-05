//! FASTA ingestion pipeline.
//!
//! Reads one or more FASTA records, 2-bit-encodes each sequence, and stores
//! the result in a [`SeqStore`].

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use scies_bio_th::io::FastaReader;
use scies_bio_th::seq::compute_sha256_hex;

use crate::model::genome::{Assembly, Chromosome, SequenceRef};
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
    /// Computes SHA-256 for each sequence and persists records to both `sequences` and `assemblies`
    /// SQLite tables when `sql_store` is supplied, along with legacy `chromosomes` table entries.
    /// Returns the list of sequence identifiers that were successfully ingested.
    pub fn ingest_file(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
        sql_store: Option<&SqliteStore>,
    ) -> Result<Vec<String>> {
        let refs = self.ingest_file_with_provenance(
            path,
            tax_id,
            assembly,
            Some(&path.to_string_lossy()),
            sql_store,
        )?;
        Ok(refs.into_iter().map(|r| r.seq_id).collect())
    }

    /// Parses `path` as a FASTA file with an explicit `source_uri`, storing 2-bit encoded sequences
    /// and recording cryptographic SHA-256 provenance in SQLite.
    pub fn ingest_file_with_provenance(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
        source_uri: Option<&str>,
        sql_store: Option<&SqliteStore>,
    ) -> Result<Vec<SequenceRef>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let fasta_reader = FastaReader::new(reader)?;

        let mut ingested_refs = Vec::new();

        if let Some(store) = sql_store {
            let asm = Assembly {
                tax_id,
                assembly_name: assembly.to_string(),
                source: source_uri.map(|s| s.to_string()).or_else(|| Some("FASTA Ingestion".to_string())),
                created_at: None,
                notes: None,
            };
            store.insert_assembly(&asm)?;
        }

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
            let sha256 = compute_sha256_hex(&record.sequence);

            self.seq_store
                .store(tax_id, assembly, &chr_name, &record.sequence)?;

            let relative_path = format!("seq/{tax_id}/{assembly}/{chr_name}.seq");
            let seq_ref = SequenceRef {
                tax_id,
                assembly_name: assembly.to_string(),
                seq_id: chr_name.clone(),
                length_bp: seq_len,
                path: relative_path,
                sha256,
                source_uri: source_uri.map(|s| s.to_string()),
            };

            if let Some(store) = sql_store {
                store.insert_sequence(&seq_ref)?;

                let is_mito = chr_name.eq_ignore_ascii_case("chrM")
                    || chr_name.eq_ignore_ascii_case("MT")
                    || chr_name.eq_ignore_ascii_case("chrMT");

                let chrom = Chromosome {
                    name: chr_name.clone(),
                    length: seq_len,
                    assembly_name: assembly.to_string(),
                    is_mitochondrial: is_mito,
                };
                let _ = store.insert_chromosome(&chrom);
            }

            ingested_refs.push(seq_ref);
        }

        Ok(ingested_refs)
    }

    /// Ingests multiple FASTA files sequentially, or concurrently when the `parallel` feature is enabled.
    ///
    /// In parallel mode, Rayon is used to distribute encoding and disk I/O across CPU cores.
    pub fn ingest_files<P: AsRef<Path> + Sync>(
        &self,
        paths: &[P],
        tax_id: u32,
        assembly: &str,
        sql_store: Option<&SqliteStore>,
    ) -> Result<Vec<String>> {
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            let results: Result<Vec<Vec<String>>> = paths
                .par_iter()
                .map(|p| self.ingest_file(p.as_ref(), tax_id, assembly, None))
                .collect();
            let all_chrs = results?.into_iter().flatten().collect::<Vec<_>>();

            if let Some(store) = sql_store {
                for chr_name in &all_chrs {
                    if let Ok(seq_len) = self.seq_store.sequence_length(tax_id, assembly, chr_name)
                    {
                        let is_mito = chr_name.eq_ignore_ascii_case("chrM")
                            || chr_name.eq_ignore_ascii_case("MT")
                            || chr_name.eq_ignore_ascii_case("chrMT");
                        let chrom = Chromosome {
                            name: chr_name.clone(),
                            length: seq_len,
                            assembly_name: assembly.to_string(),
                            is_mitochondrial: is_mito,
                        };
                        let _ = store.insert_chromosome(&chrom);
                    }
                }
            }
            Ok(all_chrs)
        }

        #[cfg(not(feature = "parallel"))]
        {
            let mut all_chrs = Vec::new();
            for p in paths {
                let chrs = self.ingest_file(p.as_ref(), tax_id, assembly, sql_store)?;
                all_chrs.extend(chrs);
            }
            Ok(all_chrs)
        }
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

    #[test]
    fn test_fasta_ingestion_with_provenance_and_sqlite() {
        let dir = tempdir().unwrap();
        let seq_dir = dir.path().join("seq");
        let meta_db = dir.path().join("meta.db");

        let seq_store = Arc::new(SeqStore::new(&seq_dir).unwrap());
        let sql_store = SqliteStore::open(&meta_db).unwrap();
        sql_store.init_schema().unwrap();

        let ingestor = FastaIngestor::new(seq_store.clone());

        let mut fasta_file = NamedTempFile::new().unwrap();
        writeln!(fasta_file, ">chr5 T2T-CHM13v2.0 TERT locus").unwrap();
        writeln!(fasta_file, "CCGTGCGTCCCAGGGCACGCACACCAGGCA").unwrap();
        fasta_file.flush().unwrap();

        let seq_refs = ingestor
            .ingest_file_with_provenance(
                fasta_file.path(),
                9606,
                "T2T-CHM13v2.0",
                Some("https://example.org/t2t.fa"),
                Some(&sql_store),
            )
            .unwrap();

        assert_eq!(seq_refs.len(), 1);
        let sref = &seq_refs[0];
        assert_eq!(sref.seq_id, "chr5");
        assert_eq!(sref.length_bp, 30);
        assert_eq!(sref.source_uri.as_deref(), Some("https://example.org/t2t.fa"));
        assert_eq!(
            sref.sha256,
            compute_sha256_hex(b"CCGTGCGTCCCAGGGCACGCACACCAGGCA")
        );

        // Verify stored in SQLite sequences table
        let from_db = sql_store
            .get_sequence_ref(9606, "T2T-CHM13v2.0", "chr5")
            .unwrap()
            .expect("sequence ref must exist in sqlite");
        assert_eq!(from_db.sha256, sref.sha256);
        assert_eq!(from_db.length_bp, 30);

        // Verify stored in SQLite assemblies table
        let asm = sql_store
            .get_assembly(9606, "T2T-CHM13v2.0")
            .unwrap()
            .expect("assembly must exist");
        assert_eq!(asm.assembly_name, "T2T-CHM13v2.0");
    }

    #[test]
    fn test_ingest_testdata_fixture_files() {
        let dir = tempdir().unwrap();
        let seq_dir = dir.path().join("seq");
        let meta_db = dir.path().join("meta.db");

        let seq_store = Arc::new(SeqStore::new(&seq_dir).unwrap());
        let sql_store = SqliteStore::open(&meta_db).unwrap();
        sql_store.init_schema().unwrap();

        let ingestor = FastaIngestor::new(seq_store.clone());

        // Test with testdata/mini_chr5_tert.fasta
        let tert_path = Path::new("testdata/mini_chr5_tert.fasta");
        if tert_path.exists() {
            let chrs = ingestor
                .ingest_file(tert_path, 9606, "T2T-CHM13v2.0", Some(&sql_store))
                .unwrap();
            assert_eq!(chrs, vec!["chr5"]);

            let sref = sql_store
                .get_sequence_ref(9606, "T2T-CHM13v2.0", "chr5")
                .unwrap()
                .unwrap();
            assert!(sref.length_bp > 0);

            let fetched = seq_store
                .fetch(9606, "T2T-CHM13v2.0", "chr5", 0, 10)
                .unwrap();
            assert_eq!(&fetched, b"CCGTGCGTCC");
        }
    }
}
