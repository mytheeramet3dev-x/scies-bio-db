//! Data model types for the biological database engine.
//!
//! All structs and enums here are serialisable and represent domain concepts
//! shared across storage, query, and API layers.

pub mod gene;
pub mod genome;
pub mod protein;
pub mod taxonomy;
pub mod variant;

// Commonly-used re-exports for ergonomic imports.
pub use gene::{Biotype, Exon, Gene, Transcript};
pub use genome::{AssemblyStatus, Chromosome, ReferenceGenome, Strand};
pub use protein::{Domain, GoAnnotation, ProteinEntry, ProteinIsoform, Ptm};
pub use taxonomy::{Lineage, Organism, TaxRank};
pub use variant::{ClinicalVariant, EffectType, VariantEffect, VariantRecord};
