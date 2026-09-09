//! Reference genome and chromosome types.

use serde::{Deserialize, Serialize};

/// Assembly status according to NCBI / INSDC classifications.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AssemblyStatus {
    /// Primary reference assembly.
    Primary,
    /// Alternate haplotype assembly.
    Alternate,
    /// Patch / fix release.
    Patch,
}

/// Strand orientation of a genomic feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Strand {
    /// Watson (plus, 5′→3′) strand.
    Forward,
    /// Crick (minus, 3′→5′) strand.
    Reverse,
}

/// Top-level description of a genome assembly.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceGenome {
    /// Assembly name (e.g. "GRCh38").
    pub assembly_name: String,
    /// NCBI Taxonomy ID of the source organism.
    pub tax_id: u32,
    /// NCBI assembly status.
    pub status: AssemblyStatus,
    /// Year the assembly was first released.
    pub release_year: u16,
}

/// A single chromosome (or scaffold) within a reference assembly.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Chromosome {
    /// Chromosome name (e.g. "chr1", "chrX", "chrM").
    pub name: String,
    /// Length in base-pairs.
    pub length: u64,
    /// Parent assembly name.
    pub assembly_name: String,
    /// `true` for the mitochondrial chromosome.
    pub is_mitochondrial: bool,
}

impl Chromosome {
    /// Returns `true` when this chromosome is an autosome (numeric name).
    ///
    /// Chromosomes whose names, after stripping an optional "chr" prefix,
    /// parse as an integer are considered autosomes.
    pub fn is_autosome(&self) -> bool {
        let name = self.name.strip_prefix("chr").unwrap_or(&self.name);
        name.parse::<u32>().is_ok()
    }

    /// Returns `true` when this chromosome is a sex chromosome (X or Y).
    pub fn is_sex_chromosome(&self) -> bool {
        let name = self.name.strip_prefix("chr").unwrap_or(&self.name);
        matches!(name, "X" | "Y" | "x" | "y")
    }
}
