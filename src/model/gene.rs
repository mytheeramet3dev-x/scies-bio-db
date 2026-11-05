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

/// A standardized gene / genomic feature record (Section 4.1).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeneFeature {
    /// NCBI Taxonomy ID of the host organism.
    pub tax_id: u32,
    /// Reference genome assembly name (e.g. "T2T-CHM13v2.0" or "GRCh38").
    pub assembly_name: String,
    /// Stable feature or gene identifier (e.g. "ENSG00000164362", "gene-TERT").
    pub gene_id: String,
    /// Gene symbol or name (e.g. "TERT", "SIRT6", "TP53").
    pub symbol: String,
    /// Chromosome / sequence name (e.g. "chr5", "chr17").
    pub chrom: String,
    /// 0-based inclusive start coordinate in base-pairs.
    pub start_bp: u64,
    /// 0-based exclusive end coordinate in base-pairs.
    pub end_bp: u64,
    /// Strand orientation (Forward or Reverse).
    pub strand: super::genome::Strand,
    /// Type of genomic feature ("gene", "transcript", "exon", etc.).
    pub feature_type: String,
    /// Data source or provider (e.g. "Ensembl", "RefSeq", "GFF3").
    pub source: Option<String>,
}

impl From<Gene> for GeneFeature {
    fn from(g: Gene) -> Self {
        Self {
            tax_id: g.tax_id,
            assembly_name: g.assembly_name,
            gene_id: g.gene_id,
            symbol: g.symbol,
            chrom: g.chr,
            start_bp: g.start,
            end_bp: g.end,
            strand: g.strand,
            feature_type: "gene".to_string(),
            source: None,
        }
    }
}
