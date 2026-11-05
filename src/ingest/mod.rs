//! Data ingestion pipelines for FASTA, GFF3, and VCF formats.

pub mod fasta;
pub mod gff3;
pub mod package;
pub mod vcf;

pub use fasta::FastaIngestor;
pub use gff3::{Gff3IngestReport, Gff3Ingestor};
pub use package::{BulkLoadReport, BulkLoader, GenomePackage};
pub use vcf::VcfIngestor;
