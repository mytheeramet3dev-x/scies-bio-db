//! Regulatory network records, multi-omics locus tracks, and spatial transcriptomics models
//! bridging `scies-bio-db` (v0.4.0) with `scies-bio-th` (v0.8.0).

use scies_bio_th::regulatory::{CorePromoterElement, RegulationType, TfbsHit};
use serde::{Deserialize, Serialize};

/// Directed Transcription Factor -> Target Gene regulatory interaction stored in SQLite.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RegulatoryEdgeRecord {
    /// Unique edge identifier (e.g. `"HUMAN_OCT4_SOX2"`).
    pub edge_id: String,
    /// NCBI Taxonomy ID (e.g. `9606` for Homo sapiens).
    pub tax_id: u32,
    /// Upstream transcription factor gene symbol.
    pub tf_symbol: String,
    /// Downstream regulated target gene symbol.
    pub target_symbol: String,
    /// Activation or repression mode.
    pub interaction: RegulationType,
    /// Regulatory coupling strength in $(0.0, 1.0]$.
    pub weight: f64,
    /// Experimental or database source (e.g. `"JASPAR-2026"`, `"ENCODE-ChIP"`).
    pub evidence_source: String,
}

/// Multi-omics track record at a gene locus stored in SQLite (`multiomics_tracks`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MultiOmicsTrackRecord {
    /// NCBI Taxonomy ID.
    pub tax_id: u32,
    /// Reference assembly name (e.g. `"GRCh38"`).
    pub assembly_name: String,
    /// Biosample or donor cohort identifier (e.g. `"Fibroblast_Donor_62yo"`).
    pub sample_id: String,
    /// HGNC gene symbol.
    pub gene_symbol: String,
    /// RNA-seq expression in Transcripts Per Million ($\text{TPM}$).
    pub rna_tpm: f64,
    /// ATAC-seq open chromatin accessibility fraction in $[0.0, 1.0]$.
    pub atac_accessibility: f64,
    /// Promoter CpG island methylation $\beta$-value in $[0.0, 1.0]$.
    pub promoter_methylation_beta: f64,
    /// H3K27ac ChIP-seq signal at the proximal enhancer.
    pub h3k27ac_chip: f64,
    /// H3K27ac input control signal.
    pub h3k27ac_input: f64,
    /// Distance in base pairs from the enhancer peak to the target TSS.
    pub enhancer_distance_bp: u64,
}

/// Spatial transcriptomics segmented cell record stored in SQLite (`spatial_cells`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpatialCellRecord {
    /// Spatial experiment dataset ID (e.g. `"MERFISH_Human_Cortex_01"`).
    pub dataset_id: String,
    /// Unique cell identifier within the section.
    pub cell_id: String,
    /// X coordinate in micrometers ($\mu\text{m}$).
    pub x_um: f64,
    /// Y coordinate in micrometers ($\mu\text{m}$).
    pub y_um: f64,
}

/// Spatial transcriptomics mRNA molecule spot stored in SQLite (`spatial_transcripts`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpatialTranscriptRecord {
    /// Spatial experiment dataset ID.
    pub dataset_id: String,
    /// Unique spot index or barcode ID.
    pub spot_id: String,
    /// Decoded gene symbol.
    pub gene_symbol: String,
    /// X coordinate in micrometers ($\mu\text{m}$).
    pub x_um: f64,
    /// Y coordinate in micrometers ($\mu\text{m}$).
    pub y_um: f64,
}

/// Combined cis-regulatory scan report over a stored promoter sequence.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PromoterScanReport {
    /// Core Pol II promoter elements detected (`TATA-box`, `Inr`, `BREu`, `DPE`).
    pub core_elements: Vec<CorePromoterElement>,
    /// JASPAR transcription factor binding site hits (`TFAP2A`, `SP1`, etc.).
    pub tfbs_hits: Vec<TfbsHit>,
}

/// Single-cell transcriptomic state record stored in SQLite (`single_cell_atlas`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SingleCellAtlasRecord {
    /// Single-cell atlas dataset identifier (e.g. `"Human_Reprogramming_Atlas_01"`).
    pub dataset_id: String,
    /// Single-cell barcode or ID.
    pub cell_id: String,
    /// Annotated cell-type or state label.
    pub cell_type: String,
    /// Expression or embedding vector across marker genes.
    pub expression: Vec<f64>,
}

/// Pharmacokinetic & pharmacodynamic compound record stored in SQLite (`pharma_compounds`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PharmaCompoundRecord {
    /// Unique compound identifier (e.g. `"CMPD_FISETIN"`, `"CMPD_NAVITOCLAX"`).
    pub compound_id: String,
    /// Compound name.
    pub compound_name: String,
    /// Primary molecular target gene symbol (e.g. `"BCL2L1"`).
    pub target_gene: String,
    /// Half-maximal effective concentration $\text{EC}_{50}$ ($\mu\text{M}$).
    pub ec50_um: f64,
    /// Median toxic concentration $\text{TD}_{50}$ ($\mu\text{M}$).
    pub td50_um: f64,
    /// Hill cooperativity coefficient $n_H$.
    pub hill_coefficient: f64,
    /// Elimination half-life $t_{1/2}$ in hours.
    pub half_life_hours: f64,
    /// Bioavailability fraction $F \in (0.0, 1.0]$.
    pub bioavailability: f64,
}

/// Synthetic gamete (MII Oocyte or Spermatozoon) generated via `rOvary` IVG and stored in SQLite (`synthetic_gametes`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SyntheticGameteRecord {
    /// Unique gamete construct identifier (e.g. `"OOCYTE_HUMAN_ROVARY_01"`).
    pub gamete_id: String,
    /// NCBI Taxonomy ID of the species (`9606` for Human).
    pub tax_id: u32,
    /// Gamete type (`"Oocyte"` or `"Spermatozoon"`).
    pub gamete_type: String,
    /// Final follicle stage achieved (e.g. `"PreOvulatoryMII"`).
    pub follicle_stage: String,
    /// Oocyte diameter in micrometers ($\mu\text{m}$).
    pub diameter_um: f64,
    /// Cytoplasmic mtDNA copy number.
    pub mitochondrial_count: u64,
    /// Maternal/paternal ICR imprint fidelity in $[0.0, 1.0]$.
    pub imprint_fidelity: f64,
    /// Estimated chromosomal aneuploidy risk in $[0.0, 1.0]$.
    pub aneuploidy_risk: f64,
    /// Overall composite gamete quality score in $[0.0, 100.0]$.
    pub overall_score: f64,
}
