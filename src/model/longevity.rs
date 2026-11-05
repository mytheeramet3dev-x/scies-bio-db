//! Longevity target catalog and genomic locus models.
//!
//! Provides models for tracking longevity-associated genes, intervention priorities,
//! and extracted genomic loci with cryptographic provenance for CRISPR/base-editing design.

use serde::{Deserialize, Serialize};
use super::gene::GeneFeature;

/// A curated longevity intervention target gene.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LongevityTarget {
    /// Standard gene symbol (e.g. "TERT", "SIRT6", "TP53").
    pub symbol: String,
    /// NCBI Taxonomy ID of the organism (e.g. 9606 for Homo sapiens).
    pub tax_id: u32,
    /// Intervention priority (1 = highest priority primary intervention target, 2 = secondary, etc.).
    pub priority: u32,
    /// Functional category (e.g. "telomere", "dna_repair", "metabolic", "tumor_suppressor", "proteostasis").
    pub category: String,
    /// Optional default or preferred assembly name (e.g. "T2T-CHM13v2.0" or "GRCh38").
    pub default_assembly: Option<String>,
    /// Optional clinical / experimental notes or biological mechanism.
    pub notes: Option<String>,
}

/// A fully-resolved genomic locus for a longevity target, combining gene annotation,
/// exact nucleotide sequence extracted from 2-bit store, and cryptographic provenance hash.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LongevityLocus {
    /// Gene annotation feature for this locus.
    pub gene: GeneFeature,
    /// Target metadata if registered in the longevity catalog.
    pub target: Option<LongevityTarget>,
    /// Exact DNA sequence extracted from the 2-bit store.
    pub sequence: scies_bio_th::Dna,
    /// 0-based coordinate start of the sequence window (inclusive).
    pub window_start: u64,
    /// 0-based coordinate end of the sequence window (exclusive).
    pub window_end: u64,
    /// Cryptographic SHA-256 hash of the parent chromosome/contig for audit provenance.
    pub provenance_sha256: String,
    /// NCBI Taxonomy ID of the host organism.
    pub tax_id: u32,
    /// Reference genome assembly name.
    pub assembly: String,
}

impl LongevityLocus {
    /// Returns the length in base-pairs of the extracted locus sequence.
    pub fn len(&self) -> usize {
        self.sequence.as_bytes().len()
    }

    /// Returns `true` if the locus sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.sequence.as_bytes().is_empty()
    }
}
