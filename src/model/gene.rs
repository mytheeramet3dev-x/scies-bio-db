//! Gene, transcript, exon, CDS, and UTR types.

use serde::{Deserialize, Serialize};

/// Biotype classification for a gene locus.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Biotype {
    /// Protein-coding gene.
    ProteinCoding,
    /// Long non-coding RNA gene.
    LongNonCodingRna,
    /// MicroRNA gene.
    MicroRna,
    /// Small nuclear RNA gene.
    SmallNuclearRna,
    /// Pseudogene locus.
    Pseudogene,
    /// Ribosomal RNA gene.
    RibosomalRna,
    /// Any other biotype not captured above.
    Other(String),
}

/// A gene locus record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    /// Stable gene identifier (Ensembl ENSG… or HGNC ID).
    pub gene_id: String,
    /// Reference genome assembly name (e.g. "GRCh38").
    pub assembly_name: String,
    /// HGNC-approved gene symbol (e.g. "TP53").
    pub symbol: String,
    /// Full gene name.
    pub name: String,
    /// Biotype of this locus.
    pub biotype: Biotype,
    /// Chromosome / contig name.
    pub chr: String,
    /// 0-based start coordinate.
    pub start: u64,
    /// 0-based exclusive end coordinate.
    pub end: u64,
    /// Strand orientation.
    pub strand: super::genome::Strand,
    /// NCBI Taxonomy ID of the host organism.
    pub tax_id: u32,
}

/// Biotype classification for an individual transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TranscriptBiotype {
    /// Messenger RNA — protein-coding transcript.
    MRNA,
    /// Non-coding RNA transcript.
    NCRNA,
    /// Pseudogene transcript.
    Pseudogene,
    /// Any other biotype.
    Other(String),
}

/// A processed transcript (mRNA or ncRNA) derived from a gene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Transcript {
    /// Stable transcript identifier (Ensembl ENST…).
    pub transcript_id: String,
    /// Parent gene identifier.
    pub gene_id: String,
    /// Reference genome assembly name.
    pub assembly_name: String,
    /// Biotype of this transcript.
    pub biotype: TranscriptBiotype,
    /// Chromosome / contig name.
    pub chr: String,
    /// 0-based start coordinate.
    pub start: u64,
    /// 0-based exclusive end coordinate.
    pub end: u64,
    /// Strand orientation.
    pub strand: super::genome::Strand,
    /// Total number of exons in this transcript.
    pub exon_count: u32,
}

/// A single exon interval belonging to a transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Exon {
    /// Stable exon identifier.
    pub exon_id: String,
    /// Parent transcript identifier.
    pub transcript_id: String,
    /// Reference genome assembly name.
    pub assembly_name: String,
    /// Chromosome / contig name.
    pub chr: String,
    /// 0-based start coordinate.
    pub start: u64,
    /// 0-based exclusive end coordinate.
    pub end: u64,
    /// Strand orientation.
    pub strand: super::genome::Strand,
    /// Position of this exon within its transcript (1-based).
    pub exon_number: u32,
}

/// A coding sequence (CDS) segment of a transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdsRegion {
    /// Parent transcript identifier.
    pub transcript_id: String,
    /// Chromosome / contig name.
    pub chr: String,
    /// 0-based start coordinate.
    pub start: u64,
    /// 0-based exclusive end coordinate.
    pub end: u64,
    /// Reading frame phase at the start of this CDS segment (0, 1, or 2).
    pub phase: u8,
}

/// An untranslated region (UTR) segment of a transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UtrRegion {
    /// Parent transcript identifier.
    pub transcript_id: String,
    /// Chromosome / contig name.
    pub chr: String,
    /// 0-based start coordinate.
    pub start: u64,
    /// 0-based exclusive end coordinate.
    pub end: u64,
    /// `true` for 5′-UTR; `false` for 3′-UTR.
    pub is_five_prime: bool,
}
