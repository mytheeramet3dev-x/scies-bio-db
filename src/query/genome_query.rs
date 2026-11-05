//! Genome-level queries: sequence fetch and chromosome listing.

use std::sync::Arc;

use crate::index::KmerIndex;
use crate::model::genome::{Chromosome, Strand};
use crate::storage::{SeqStore, SqliteStore};
use crate::{BioDbError, Result};

/// A located CRISPR off-target site discovered across a chromosome in `BioDb`.
#[derive(Debug, Clone, PartialEq)]
pub struct CrisprOffTargetHit {
    /// Chromosome identifier.
    pub chr: String,
    /// 0-based start coordinate of the 20-nt protospacer on the forward strand.
    pub position: u64,
    /// Strand where the PAM-adjacent match occurs.
    pub strand: Strand,
    /// 3-nt PAM sequence (e.g. `AGG`, `TGG`, `NAG`).
    pub pam: String,
    /// Number of mismatches against the query 20-nt guide RNA.
    pub mismatches: usize,
    /// Doench 2016 Cutting Frequency Determination (CFD) cleavage score (0.0 to 100.0).
    pub cfd_score: f64,
}

/// Multi-chromosome / assembly-wide CRISPR SpCas9 off-target scan report.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyOffTargetReport {
    /// Assembly name (e.g. `"GRCh38.p14"` or `"T2T-CHM13v2.0"`).
    pub assembly: String,
    /// Number of chromosomes scanned in the database for this assembly.
    pub chromosomes_scanned: usize,
    /// Total base pairs scanned across all stored chromosomes in the assembly.
    pub total_bp_scanned: u64,
    /// Count of exact 0-mismatch on-target sites (`NGG` / `NAG`).
    pub exact_matches: usize,
    /// All discovered candidate cleavage sites sorted by descending `cfd_score`.
    pub hits: Vec<CrisprOffTargetHit>,
    /// Cumulative CFD cleavage burden across all off-target sites (`mismatches > 0` or secondary hits).
    pub aggregate_offtarget_cfd_burden: f64,
    /// Aggregate specificity score ($0.0$ to $100.0$), defined as $\frac{100}{100 + \sum \text{CFD}_{\text{off}}} \times 100$.
    pub specificity_score: f64,
}

/// Handles queries against chromosome sequences and assembly metadata.
pub struct GenomeQuery {
    /// Shared sequence storage backend.
    pub(crate) seq_store: Arc<SeqStore>,
    /// Shared SQLite metadata store.
    pub(crate) sql_store: Arc<SqliteStore>,
}

impl GenomeQuery {
    /// Creates a new [`GenomeQuery`].
    pub fn new(seq_store: Arc<SeqStore>, sql_store: Arc<SqliteStore>) -> Self {
        Self {
            seq_store,
            sql_store,
        }
    }

    /// Fetches the raw DNA sequence for `chr` over the half-open interval
    /// `[start, end)` and returns it as a [`scies_bio_th::Dna`] value.
    pub fn fetch_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<scies_bio_th::Dna> {
        let raw_bytes = self.seq_store.fetch(tax_id, assembly, chr, start, end)?;
        scies_bio_th::Dna::new(&raw_bytes).map_err(BioDbError::Bio)
    }

    /// Fetches raw sequence bytes for `chr` over the half-open interval `[start, end)`.
    ///
    /// Preserves full IUPAC degenerate symbols and Uracil without converting to [`scies_bio_th::Dna`].
    pub fn fetch_raw_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>> {
        self.seq_store.fetch(tax_id, assembly, chr, start, end)
    }

    /// Fetches the full DNA sequence for `chr` and returns it as a [`scies_bio_th::Dna`] value.
    pub fn fetch_full_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
    ) -> Result<scies_bio_th::Dna> {
        let len = self.seq_store.sequence_length(tax_id, assembly, chr)?;
        self.fetch_sequence(tax_id, assembly, chr, 0, len)
    }

    /// Fetches the complete raw sequence bytes for `chr`.
    pub fn fetch_full_raw_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
    ) -> Result<Vec<u8>> {
        let len = self.seq_store.sequence_length(tax_id, assembly, chr)?;
        self.fetch_raw_sequence(tax_id, assembly, chr, 0, len)
    }

    /// Returns metadata for every chromosome in assembly `assembly`.
    pub fn list_chromosomes(&self, assembly: &str) -> Result<Vec<Chromosome>> {
        self.sql_store.list_chromosomes(assembly)
    }

    /// Builds an in-memory [`KmerIndex`] of length `k` for chromosome `chr`.
    pub fn build_kmer_index(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        k: usize,
    ) -> Result<KmerIndex> {
        let raw = self.fetch_full_raw_sequence(tax_id, assembly, chr)?;
        let mut index = KmerIndex::new(k);
        index.index_sequence(chr, &raw);
        Ok(index)
    }

    /// Builds an in-memory [`crate::index::FmIndex`] for chromosome `chr` supporting
    /// $O(|P|)$ exact backward search and bounded-mismatch approximate alignment.
    pub fn build_fm_index(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
    ) -> Result<crate::index::FmIndex> {
        let raw = self.fetch_full_raw_sequence(tax_id, assembly, chr)?;
        Ok(crate::index::FmIndex::build(chr, &raw))
    }

    /// Builds an in-memory assembly-wide [`crate::index::MultiChromFmIndex`] covering all chromosomes
    /// stored under `assembly`.
    pub fn build_assembly_fm_index(
        &self,
        tax_id: u32,
        assembly: &str,
    ) -> Result<crate::index::MultiChromFmIndex> {
        let chroms = self.list_chromosomes(assembly)?;
        let mut multi = crate::index::MultiChromFmIndex::new();
        for chrom in &chroms {
            let raw = self.fetch_full_raw_sequence(tax_id, assembly, &chrom.name)?;
            multi.add_chromosome(&chrom.name, &raw);
        }
        Ok(multi)
    }

    /// Computes the Jaccard k-mer similarity between chromosome `chr` and `query_sequence`.
    pub fn estimate_similarity(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        query_sequence: &[u8],
        k: usize,
    ) -> Result<f64> {
        let index = self.build_kmer_index(tax_id, assembly, chr, k)?;
        Ok(index.jaccard_similarity(query_sequence))
    }

    /// Scans a chromosome in `BioDb` for potential SpCas9 (`NGG` / `NAG`) off-target sites
    /// within `max_mismatches` of `guide_rna`, scoring each hit with the Doench 2016 CFD matrix.
    pub fn scan_crispr_offtargets(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        guide_rna: &scies_bio_th::Rna,
        max_mismatches: usize,
    ) -> Result<Vec<CrisprOffTargetHit>> {
        let dna = self.fetch_full_sequence(tax_id, assembly, chr)?;
        let mut hits = Vec::new();
        let guide_len = guide_rna.sequence.len();
        if guide_len == 0 || dna.sequence.len() < guide_len + 3 {
            return Ok(hits);
        }

        // 1. Scan Forward strand
        let fwd_bytes = &dna.sequence;
        for pam_idx in guide_len..=(fwd_bytes.len() - 3) {
            let p1 = fwd_bytes[pam_idx + 1];
            let p2 = fwd_bytes[pam_idx + 2];
            if (p1 == b'G' || p1 == b'A') && p2 == b'G' {
                let start = pam_idx - guide_len;
                if let Ok(candidate_dna) = scies_bio_th::Dna::new(&fwd_bytes[start..pam_idx]) {
                    let cfd = scies_bio_th::edit::offtarget::calculate_doench_2016_cfd(
                        guide_rna,
                        &candidate_dna,
                    );
                    if cfd.mismatches <= max_mismatches {
                        let pam_weight = if p1 == b'G' { 1.0 } else { 0.25 };
                        hits.push(CrisprOffTargetHit {
                            chr: chr.to_string(),
                            position: start as u64,
                            strand: Strand::Forward,
                            pam: String::from_utf8_lossy(&fwd_bytes[pam_idx..pam_idx + 3])
                                .to_string(),
                            mismatches: cfd.mismatches,
                            cfd_score: cfd.score * pam_weight,
                        });
                    }
                }
            }
        }

        // 2. Scan Reverse strand
        let rev_dna = dna.reverse_complement();
        let rev_bytes = &rev_dna.sequence;
        let total_len = fwd_bytes.len();
        for pam_idx in guide_len..=(rev_bytes.len() - 3) {
            let p1 = rev_bytes[pam_idx + 1];
            let p2 = rev_bytes[pam_idx + 2];
            if (p1 == b'G' || p1 == b'A') && p2 == b'G' {
                let start = pam_idx - guide_len;
                if let Ok(candidate_dna) = scies_bio_th::Dna::new(&rev_bytes[start..pam_idx]) {
                    let cfd = scies_bio_th::edit::offtarget::calculate_doench_2016_cfd(
                        guide_rna,
                        &candidate_dna,
                    );
                    if cfd.mismatches <= max_mismatches {
                        let pam_weight = if p1 == b'G' { 1.0 } else { 0.25 };
                        let fwd_pos = total_len.saturating_sub(pam_idx) as u64;
                        hits.push(CrisprOffTargetHit {
                            chr: chr.to_string(),
                            position: fwd_pos,
                            strand: Strand::Reverse,
                            pam: String::from_utf8_lossy(&rev_bytes[pam_idx..pam_idx + 3])
                                .to_string(),
                            mismatches: cfd.mismatches,
                            cfd_score: cfd.score * pam_weight,
                        });
                    }
                }
            }
        }

        hits.sort_by(|a, b| {
            b.cfd_score
                .partial_cmp(&a.cfd_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Ok(hits)
    }

    /// Scans all chromosomes registered under `assembly` in `BioDb` for SpCas9 (`NGG` / `NAG`)
    /// on-target and off-target sites within `max_mismatches` of `guide_rna`, computing an
    /// assembly-wide specificity score and off-target CFD burden.
    pub fn scan_assembly_crispr_offtargets(
        &self,
        tax_id: u32,
        assembly: &str,
        guide_rna: &scies_bio_th::Rna,
        max_mismatches: usize,
    ) -> Result<AssemblyOffTargetReport> {
        let chroms = self.list_chromosomes(assembly)?;
        let mut all_hits = Vec::new();
        let mut total_bp_scanned = 0u64;

        for chrom in &chroms {
            total_bp_scanned += chrom.length;
            let chr_hits = self.scan_crispr_offtargets(
                tax_id,
                assembly,
                &chrom.name,
                guide_rna,
                max_mismatches,
            )?;
            all_hits.extend(chr_hits);
        }

        all_hits.sort_by(|a, b| {
            b.cfd_score
                .partial_cmp(&a.cfd_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let exact_matches = all_hits.iter().filter(|h| h.mismatches == 0).count();
        // Sum CFD across all mismatched off-targets plus any secondary exact matches beyond the first on-target site
        let mut seen_primary_exact = false;
        let mut aggregate_offtarget_cfd_burden = 0.0;
        for hit in &all_hits {
            if hit.mismatches == 0 && !seen_primary_exact {
                seen_primary_exact = true;
            } else {
                aggregate_offtarget_cfd_burden += hit.cfd_score;
            }
        }

        let specificity_score =
            (100.0 / (100.0 + aggregate_offtarget_cfd_burden) * 100.0).clamp(0.0, 100.0);

        Ok(AssemblyOffTargetReport {
            assembly: assembly.to_string(),
            chromosomes_scanned: chroms.len(),
            total_bp_scanned,
            exact_matches,
            hits: all_hits,
            aggregate_offtarget_cfd_burden,
            specificity_score,
        })
    }

    /// Scans a stored promoter region `[start, end)` on `chr` for Pol II core promoter motifs
    /// (`TATA-box`, `Inr`, `BREu`, `DPE`) and canonical JASPAR TFBS (`TFAP2A`, `SP1`).
    pub fn scan_promoter_elements(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<crate::model::PromoterScanReport> {
        let dna = self.fetch_sequence(tax_id, assembly, chr, start, end)?;
        let core_elements = scies_bio_th::regulatory::detect_core_promoter_elements(&dna, 0.85);
        let mut tfbs_hits = scies_bio_th::regulatory::scan_jaspar_tfbs(
            &dna,
            &scies_bio_th::regulatory::jaspar_ma0003_1_tfap2a(),
            0.85,
        );
        let sp1_hits = scies_bio_th::regulatory::scan_jaspar_tfbs(
            &dna,
            &scies_bio_th::regulatory::jaspar_ma0079_3_sp1(),
            0.85,
        );
        tfbs_hits.extend(sp1_hits);

        Ok(crate::model::PromoterScanReport {
            core_elements,
            tfbs_hits,
        })
    }
}
