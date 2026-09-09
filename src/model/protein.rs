//! Protein entry, isoform, domain, PTM, and GO annotation types.

use serde::{Deserialize, Serialize};

/// UniProtKB review status.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ReviewStatus {
    /// Manually reviewed (Swiss-Prot).
    SwissProt,
    /// Computationally analysed (TrEMBL / unreviewed).
    TrEMBL,
}

/// Gross subcellular location of a protein.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SubcellularLocation {
    /// Cell nucleus.
    Nucleus,
    /// Cytoplasm / cytosol.
    Cytoplasm,
    /// Plasma or intracellular membrane.
    Membrane,
    /// Mitochondrial compartment.
    Mitochondria,
    /// Endoplasmic reticulum.
    EndoplasmicReticulum,
    /// Extracellular space / secreted.
    Extracellular,
    /// Any other location.
    Other(String),
}

/// A UniProtKB protein entry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinEntry {
    /// Primary UniProtKB accession (e.g. "P04637").
    pub uniprot_accession: String,
    /// UniProtKB entry name (e.g. "P53_HUMAN").
    pub entry_name: String,
    /// HGNC gene symbol of the encoding gene.
    pub gene_symbol: String,
    /// NCBI Taxonomy ID.
    pub tax_id: u32,
    /// Review status of this entry.
    pub status: ReviewStatus,
    /// Canonical amino-acid sequence (single-letter code).
    pub sequence: String,
    /// Number of amino acids.
    pub length: u32,
    /// Known subcellular locations.
    pub subcellular_locations: Vec<SubcellularLocation>,
}

/// An alternatively-spliced protein isoform.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinIsoform {
    /// Isoform identifier (e.g. "P04637-2").
    pub isoform_id: String,
    /// Primary accession of the canonical entry.
    pub uniprot_accession: String,
    /// Amino-acid sequence of this isoform.
    pub sequence: String,
    /// Number of amino acids.
    pub length: u32,
}

/// A functional or structural domain annotated on a protein.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Domain {
    /// Database-specific domain identifier (e.g. "PF00870").
    pub domain_id: String,
    /// Human-readable domain name.
    pub name: String,
    /// Source database identifier (e.g. "pfam", "interpro").
    pub db: String,
    /// UniProtKB accession of the parent protein.
    pub uniprot_accession: String,
    /// 1-based start amino-acid position.
    pub start_aa: u32,
    /// 1-based end amino-acid position (inclusive).
    pub end_aa: u32,
}

/// Type of post-translational modification (PTM).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PtmType {
    /// Phosphorylation on Ser/Thr/Tyr.
    Phosphorylation,
    /// N- or O-glycosylation.
    Glycosylation,
    /// Ubiquitin conjugation.
    Ubiquitination,
    /// Methyl group addition.
    Methylation,
    /// Acetyl group addition.
    Acetylation,
    /// Any other PTM.
    Other(String),
}

/// A post-translational modification site annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ptm {
    /// UniProtKB accession of the modified protein.
    pub uniprot_accession: String,
    /// 1-based position of the modified residue.
    pub position: u32,
    /// Type of modification.
    pub ptm_type: PtmType,
    /// Experimental or computational evidence code.
    pub evidence: String,
}

/// Gene Ontology (GO) aspect.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GoAspect {
    /// Biological process aspect (P).
    BiologicalProcess,
    /// Molecular function aspect (F).
    MolecularFunction,
    /// Cellular component aspect (C).
    CellularComponent,
}

/// A single GO annotation linking a protein to a GO term.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoAnnotation {
    /// UniProtKB accession of the annotated protein.
    pub uniprot_accession: String,
    /// GO identifier (e.g. "GO:0005634").
    pub go_id: String,
    /// Human-readable term name.
    pub term: String,
    /// Ontology aspect of the term.
    pub aspect: GoAspect,
}
