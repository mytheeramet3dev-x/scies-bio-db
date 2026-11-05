//! Epigenomic features, 3D protein structure records, and dataset provenance manifests.

use serde::{Deserialize, Serialize};

/// Type of epigenomic or regulatory annotation track.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EpigeneticTrackType {
    /// Gardiner-Garden & Frommer CpG Island.
    CpgIsland,
    /// Histone post-translational modification mark (e.g. H3K4me3, H3K9me3, H3K27ac).
    HistoneMark,
    /// Open chromatin accessibility peak (ATAC-seq / DNase-seq).
    ChromatinPeak,
    /// Promoter or enhancer cis-regulatory element.
    RegulatoryElement,
}

/// An epigenomic or regulatory feature stored in `BioDb`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EpigeneticFeature {
    /// Unique feature identifier.
    pub feature_id: String,
    /// NCBI Taxonomy ID.
    pub tax_id: u32,
    /// Reference assembly name (e.g. "T2T-CHM13v2.0", "GRCh38").
    pub assembly_name: String,
    /// Chromosome / contig name.
    pub chr: String,
    /// 0-based start coordinate (inclusive).
    pub start: u64,
    /// 0-based end coordinate (exclusive).
    pub end: u64,
    /// Track category.
    pub track_type: EpigeneticTrackType,
    /// Label or mark name (e.g. "H3K9me3", "CpG_Island_TERT", "ATAC_Peak_1").
    pub mark_label: String,
    /// Quantitative signal or methylation beta/obs-exp score.
    pub signal_value: f64,
}

/// 3D macromolecular structure metadata linked to a protein entry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ProteinStructureRecord {
    /// PDB or AlphaFold accession identifier (e.g. "1TUP", "AF-P04637-F1").
    pub structure_id: String,
    /// Associated UniProt accession or gene symbol.
    pub gene_symbol: String,
    /// NCBI Taxonomy ID.
    pub tax_id: u32,
    /// Experimental or predicted resolution in Angstroms.
    pub resolution_angstrom: f64,
    /// Number of C-alpha residues modeled.
    pub residue_count: u32,
    /// Primary binding pocket volume in cubic Angstroms ($Å^3$).
    pub pocket_volume_a3: f64,
}

/// Cryptographic provenance and laboratory validation manifest for an ingested dataset.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DatasetManifest {
    /// Unique manifest identifier (e.g. "T2T-CHM13v2.0-2026").
    pub manifest_id: String,
    /// NCBI Taxonomy ID.
    pub tax_id: u32,
    /// Assembly accession or name.
    pub assembly_name: String,
    /// Data provider (e.g. "NCBI RefSeq", "T2T Consortium", "Ensembl").
    pub provider: String,
    /// SHA-256 hex digest of the source artifact.
    pub source_sha256: String,
    /// ISO-8601 timestamp or date of ingestion.
    pub ingested_at: String,
}
