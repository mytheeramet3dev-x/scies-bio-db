//! Multi-file genome package definition and bulk loader.
//!
//! Models standard NCBI (RefSeq/GenBank) and Ensembl release packages containing
//! FASTA genomic sequences, GFF3 gene annotations, and VCF variant sets.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::{FastaIngestor, Gff3Ingestor, VcfIngestor};
use crate::model::genome::ReferenceGenome;
use crate::model::taxonomy::Organism;
use crate::storage::{SeqStore, SqliteStore};
use crate::Result;

/// Complete metadata and file paths representing an NCBI or Ensembl biological dataset release.
#[derive(Debug, Clone)]
pub struct GenomePackage {
    /// Organism taxonomy metadata.
    pub organism: Organism,
    /// Reference genome assembly metadata.
    pub reference_genome: ReferenceGenome,
    /// One or more FASTA sequence file paths (e.g. chromosomes, contigs, or whole genome).
    pub fasta_files: Vec<PathBuf>,
    /// Optional GFF3 feature annotation file paths.
    pub gff3_files: Vec<PathBuf>,
    /// Optional VCF variant call file paths.
    pub vcf_files: Vec<PathBuf>,
}

impl GenomePackage {
    /// Constructs a new [`GenomePackage`] for a given organism and reference assembly.
    pub fn new(organism: Organism, reference_genome: ReferenceGenome) -> Self {
        Self {
            organism,
            reference_genome,
            fasta_files: Vec::new(),
            gff3_files: Vec::new(),
            vcf_files: Vec::new(),
        }
    }

    /// Adds a FASTA sequence file path to the package.
    pub fn add_fasta<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.fasta_files.push(path.as_ref().to_path_buf());
        self
    }

    /// Adds a GFF3 annotation file path to the package.
    pub fn add_gff3<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.gff3_files.push(path.as_ref().to_path_buf());
        self
    }

    /// Adds a VCF variant file path to the package.
    pub fn add_vcf<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.vcf_files.push(path.as_ref().to_path_buf());
        self
    }
}

/// Cumulative statistics and summary from a bulk loading operation.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BulkLoadReport {
    /// Total chromosome records ingested.
    pub chromosomes_ingested: usize,
    /// Total nucleotide base count across all ingested chromosomes.
    pub total_bases_ingested: u64,
    /// Total genomic features (genes, transcripts, exons) ingested.
    pub features_ingested: u64,
    /// Total variant calls ingested.
    pub variants_ingested: u64,
}

/// Orchestrates atomic bulk loading of an entire [`GenomePackage`] into [`SeqStore`] and [`SqliteStore`].
pub struct BulkLoader {
    seq_store: Arc<SeqStore>,
    sql_store: Arc<SqliteStore>,
}

impl BulkLoader {
    /// Creates a new [`BulkLoader`].
    pub fn new(seq_store: Arc<SeqStore>, sql_store: Arc<SqliteStore>) -> Self {
        Self {
            seq_store,
            sql_store,
        }
    }

    /// Executes the full lifecycle ingestion of `package`:
    /// 1. Registers the organism and reference genome in the SQLite catalog.
    /// 2. Enables high-throughput bulk mode on SQLite (`PRAGMA cache_size = -64000; temp_store = MEMORY`).
    /// 3. Ingests all FASTA sequences into `.seq v2` storage (leveraging Rayon in parallel mode).
    /// 4. Ingests all GFF3 annotations in batch transactions.
    /// 5. Ingests all VCF variants in chunked batch transactions.
    /// 6. Restores standard SQLite operational settings.
    pub fn load(&self, package: &GenomePackage) -> Result<BulkLoadReport> {
        let tax_id = package.organism.tax_id;
        let assembly = &package.reference_genome.assembly_name;

        // 1. Register metadata entities
        self.sql_store.insert_organism(&package.organism)?;
        self.sql_store
            .insert_reference_genome(&package.reference_genome)?;

        let mut report = BulkLoadReport::default();

        self.sql_store.with_bulk_mode(|| {
            // 2. Ingest FASTA sequence files (in parallel if feature enabled)
            let fasta_ingestor = FastaIngestor::new(Arc::clone(&self.seq_store));
            let ingested_chrs = fasta_ingestor.ingest_files(
                &package.fasta_files,
                tax_id,
                assembly,
                Some(&self.sql_store),
            )?;
            report.chromosomes_ingested = ingested_chrs.len();
            for chr in &ingested_chrs {
                if let Ok(len) = self.seq_store.sequence_length(tax_id, assembly, chr) {
                    report.total_bases_ingested += len;
                }
            }

            // 3. Ingest GFF3 annotation files
            let gff_ingestor = Gff3Ingestor::new();
            for gff_path in &package.gff3_files {
                let count =
                    gff_ingestor.ingest_file(gff_path, tax_id, assembly, &self.sql_store)?;
                report.features_ingested += count;
            }

            // 4. Ingest VCF variant files
            let vcf_ingestor = VcfIngestor::new();
            for vcf_path in &package.vcf_files {
                let count =
                    vcf_ingestor.ingest_file(vcf_path, tax_id, assembly, &self.sql_store)?;
                report.variants_ingested += count;
            }

            Ok(())
        })?;

        Ok(report)
    }
}
