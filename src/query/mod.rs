//! Query interfaces for genome, gene, protein, and variant data.

pub mod gene_query;
pub mod genome_query;
pub mod protein_query;
pub mod variant_query;

pub use gene_query::GeneQuery;
pub use genome_query::GenomeQuery;
pub use protein_query::ProteinQuery;
pub use variant_query::VariantQuery;
