//! Genomic variant, clinical annotation, and functional effect types.

use serde::{Deserialize, Serialize};

/// Classification of the variant class by mutation type.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum VariantType {
    /// Single nucleotide polymorphism.
    Snp,
    /// Insertion of one or more bases.
    Insertion,
    /// Deletion of one or more bases.
    Deletion,
    /// Multi-nucleotide polymorphism (block substitution).
    Mnp,
    /// Large structural variant (SV/CNV).
    StructuralVariant,
}

/// ClinVar / ACMG clinical significance classification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ClinicalSignificance {
    /// Classified as disease-causing.
    Pathogenic,
    /// Likely disease-causing.
    LikelyPathogenic,
    /// Variant of uncertain significance.
    VUS,
    /// Likely neutral / benign.
    LikelyBenign,
    /// Classified as neutral / benign.
    Benign,
    /// Conflicting interpretations in the literature.
    Conflicting,
    /// Significance not determined.
    Unknown,
}

/// A raw variant record from a population or reference database.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantRecord {
    /// dbSNP rs identifier (e.g. "rs1042522").
    pub variant_id: String,
    /// Chromosome / contig name.
    pub chr: String,
    /// 1-based position on the chromosome.
    pub position: u64,
    /// Reference allele sequence.
    pub reference: String,
    /// Alternate allele sequence.
    pub alternate: String,
    /// Variant class.
    pub variant_type: VariantType,
    /// NCBI Taxonomy ID of the host organism.
    pub tax_id: u32,
}

/// Clinical interpretation of a variant from ClinVar or similar.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClinicalVariant {
    /// Variant identifier linking to [`VariantRecord`].
    pub variant_id: String,
    /// Clinical significance classification.
    pub clinical_significance: ClinicalSignificance,
    /// Associated disease / phenotype name.
    pub condition: String,
    /// ClinVar review status string.
    pub review_status: String,
    /// Affected gene symbol, if known.
    pub gene_symbol: Option<String>,
}

/// Predicted or observed functional effect of a variant on a transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum EffectType {
    /// Amino-acid changing substitution.
    Missense,
    /// Codon change preserving the amino acid.
    Synonymous,
    /// Premature stop codon introduced.
    Nonsense,
    /// Insertion / deletion causing a reading-frame shift.
    Frameshift,
    /// Variant affecting a splice site.
    Splice,
    /// Variant in the upstream promoter region.
    Upstream,
    /// Variant in the downstream region.
    Downstream,
    /// Intergenic variant (no overlapping gene).
    Intergenic,
    /// Any other effect.
    Other(String),
}

/// The functional consequence of a variant on a specific transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantEffect {
    /// Variant identifier linking to [`VariantRecord`].
    pub variant_id: String,
    /// Transcript on which the effect was computed.
    pub transcript_id: String,
    /// Parent gene of the transcript.
    pub gene_id: String,
    /// Predicted effect type.
    pub effect_type: EffectType,
    /// HGVS-style codon change string (e.g. "cgT/cgC"), if applicable.
    pub codon_change: Option<String>,
    /// HGVS-style amino-acid change string (e.g. "p.Arg72Pro"), if applicable.
    pub aa_change: Option<String>,
}
