//! Data model types for the biological database engine.
//!
//! All structs and enums here are serialisable and represent domain concepts
//! shared across storage, query, and API layers.

pub mod epigenome;
pub mod gene;
pub mod genome;
pub mod longevity;
pub mod protein;
pub mod regulatory;
pub mod taxonomy;
pub mod variant;

// Commonly-used re-exports for ergonomic imports.
pub use epigenome::{
    DatasetManifest, EpigeneticFeature, EpigeneticTrackType, ProteinStructureRecord,
};
pub use gene::{Biotype, Exon, Gene, GeneFeature, Transcript};
pub use genome::{Assembly, AssemblyStatus, Chromosome, ReferenceGenome, SequenceRef, Strand};
pub use longevity::{LongevityLocus, LongevityTarget};
pub use protein::{Domain, GoAnnotation, ProteinEntry, ProteinIsoform, Ptm};
pub use regulatory::{
    MultiOmicsTrackRecord, PharmaCompoundRecord, PromoterScanReport, RegulatoryEdgeRecord,
    SingleCellAtlasRecord, SpatialCellRecord, SpatialTranscriptRecord, SyntheticGameteRecord,
};
pub use taxonomy::{Lineage, Organism, TaxRank};
pub use variant::{ClinicalVariant, EffectType, VariantEffect, VariantRecord};
