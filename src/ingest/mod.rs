//! Data ingestion pipelines for FASTA, GFF3, and VCF formats.

pub mod fasta;
pub mod gff3;
pub mod vcf;

pub use fasta::FastaIngestor;
pub use gff3::Gff3Ingestor;
pub use vcf::VcfIngestor;
