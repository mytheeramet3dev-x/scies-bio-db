//! Burrows-Wheeler Transform (BWT) and Ferragina-Manzini (`FM-Index`) for $O(|P|)$ exact and
//! bounded-backtracking approximate sequence alignment across chromosomes.
//!
//! References:
//! - Ferragina & Manzini (2000) *41st IEEE FOCS* 390-398 (Opportunistic data structures / FM-index).
//! - Langmead et al. (2009) *Genome Biology* 10:R25 (Ultrafast memory-efficient alignment — Bowtie).
//! - Li & Durbin (2009) *Bioinformatics* 25:1754-1760 (Fast and accurate short read alignment with BWT — BWA).

use std::collections::HashMap;

/// Alphabet symbol encoding: `0 = '$'` (sentinel), `1 = 'A'`, `2 = 'C'`, `3 = 'G'`, `4 = 'T'`.
const ALPHABET_SIZE: usize = 5;

#[inline]
fn encode_base(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'$' | 0 => Some(0),
        b'A' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        b'T' | b'U' => Some(4),
        _ => None,
    }
}

/// Ferragina-Manzini (`FM-Index`) over a single chromosome sequence, supporting $O(|P|)$ exact
/// backward search and BWA-style bounded-mismatch backtracking.
#[derive(Debug, Clone)]
pub struct FmIndex {
    /// Chromosome identifier.
    pub chr: String,
    /// Original sequence length (excluding sentinel `$`).
    pub seq_len: usize,
    /// Suffix Array `SA[0..=seq_len]`.
    sa: Vec<u32>,
    /// Burrows-Wheeler Transform encoded symbols (`0..5`).
    bwt: Vec<u8>,
    /// Cumulative frequency table `C[0..5]` (count of symbols strictly smaller than `c`).
    c_table: [usize; ALPHABET_SIZE],
    /// Prefix occurrence matrix `occ[i][c]` = count of symbol `c` in `bwt[0..i]`.
    occ: Vec<[u32; ALPHABET_SIZE]>,
}

impl FmIndex {
    /// Builds an [`FmIndex`] for chromosome `chr` from raw DNA bytes.
    ///
    /// Non-ACGT IUPAC ambiguity characters are normalized to `'A'` during indexing.
    pub fn build(chr: &str, sequence: &[u8]) -> Self {
        let n = sequence.len();
        let mut text = Vec::with_capacity(n + 1);
        for &b in sequence {
            let sym = encode_base(b).unwrap_or(1) as u8;
            text.push(sym);
        }
        text.push(0u8); // Append lexicographically smallest sentinel '$' = 0

        let total = n + 1;
        let mut sa: Vec<u32> = (0..total as u32).collect();
        sa.sort_by(|&a, &b| text[a as usize..].cmp(&text[b as usize..]));

        let mut bwt = Vec::with_capacity(total);
        let mut counts = [0usize; ALPHABET_SIZE];
        let mut occ = Vec::with_capacity(total + 1);
        let mut running = [0u32; ALPHABET_SIZE];
        occ.push(running);

        for &s in &sa {
            let prev_idx = if s == 0 { total - 1 } else { (s as usize) - 1 };
            let sym = text[prev_idx];
            bwt.push(sym);
            counts[sym as usize] += 1;
            running[sym as usize] += 1;
            occ.push(running);
        }

        let mut c_table = [0usize; ALPHABET_SIZE];
        let mut cum = 0usize;
        for c in 0..ALPHABET_SIZE {
            c_table[c] = cum;
            cum += counts[c];
        }

        Self {
            chr: chr.to_string(),
            seq_len: n,
            sa,
            bwt,
            c_table,
            occ,
        }
    }

    /// Returns the raw Burrows-Wheeler Transform bytes as ASCII (`$`, `A`, `C`, `G`, `T`).
    pub fn bwt_ascii(&self) -> Vec<u8> {
        self.bwt
            .iter()
            .map(|&s| match s {
                0 => b'$',
                1 => b'A',
                2 => b'C',
                3 => b'G',
                _ => b'T',
            })
            .collect()
    }

    /// Performs $O(|P|)$ backward search and returns the Suffix Array interval `[sp, ep)`.
    pub fn backward_search_interval(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        if pattern.is_empty() {
            return None;
        }
        let mut sp = 0usize;
        let mut ep = self.sa.len();

        for &b in pattern.iter().rev() {
            let c = encode_base(b)?;
            if c == 0 {
                return None;
            }
            sp = self.c_table[c] + self.occ[sp][c] as usize;
            ep = self.c_table[c] + self.occ[ep][c] as usize;
            if sp >= ep {
                return None;
            }
        }
        Some((sp, ep))
    }

    /// Counts exact occurrences of `pattern` in $O(|P|)$ time.
    pub fn count_exact(&self, pattern: &[u8]) -> usize {
        self.backward_search_interval(pattern)
            .map(|(sp, ep)| ep - sp)
            .unwrap_or(0)
    }

    /// Locates all 0-based start coordinates of exact matches of `pattern`, sorted ascending.
    pub fn locate_exact(&self, pattern: &[u8]) -> Vec<u64> {
        let Some((sp, ep)) = self.backward_search_interval(pattern) else {
            return Vec::new();
        };
        let mut positions: Vec<u64> = self.sa[sp..ep].iter().map(|&p| p as u64).collect();
        positions.sort_unstable();
        positions
    }

    /// Locates all occurrences of `pattern` within at most `max_mismatches` substitutions using
    /// BWA-style bounded DFS backtracking over the FM-index, returning `(position, mismatches)`.
    pub fn locate_approximate(
        &self,
        pattern: &[u8],
        max_mismatches: usize,
    ) -> Vec<(u64, usize)> {
        if pattern.is_empty() || pattern.len() > self.seq_len {
            return Vec::new();
        }
        let mut pat_encoded = Vec::with_capacity(pattern.len());
        for &b in pattern {
            let Some(c) = encode_base(b) else {
                return Vec::new();
            };
            pat_encoded.push(c);
        }

        let mut best_by_pos: HashMap<u64, usize> = HashMap::new();
        self.dfs_backtrack(
            &pat_encoded,
            pattern.len(),
            0,
            self.sa.len(),
            0,
            max_mismatches,
            &mut best_by_pos,
        );

        let max_valid_start = (self.seq_len - pattern.len()) as u64;
        let mut hits: Vec<(u64, usize)> = best_by_pos
            .into_iter()
            .filter(|&(pos, _)| pos <= max_valid_start)
            .collect();
        hits.sort_by_key(|&(pos, mm)| (mm, pos));
        hits
    }

    #[allow(clippy::too_many_arguments)]
    fn dfs_backtrack(
        &self,
        pat: &[usize],
        rem_idx: usize,
        sp: usize,
        ep: usize,
        mismatches: usize,
        max_mismatches: usize,
        out: &mut HashMap<u64, usize>,
    ) {
        if sp >= ep || mismatches > max_mismatches {
            return;
        }
        if rem_idx == 0 {
            for &sa_pos in &self.sa[sp..ep] {
                let pos = sa_pos as u64;
                let entry = out.entry(pos).or_insert(mismatches);
                if mismatches < *entry {
                    *entry = mismatches;
                }
            }
            return;
        }

        let target_c = pat[rem_idx - 1];
        for c in 1..ALPHABET_SIZE {
            let next_sp = self.c_table[c] + self.occ[sp][c] as usize;
            let next_ep = self.c_table[c] + self.occ[ep][c] as usize;
            if next_sp < next_ep {
                let next_mm = mismatches + usize::from(c != target_c);
                if next_mm <= max_mismatches {
                    self.dfs_backtrack(
                        pat,
                        rem_idx - 1,
                        next_sp,
                        next_ep,
                        next_mm,
                        max_mismatches,
                        out,
                    );
                }
            }
        }
    }
}

/// Multi-chromosome FM-Index collection covering an entire genome assembly.
#[derive(Debug, Clone, Default)]
pub struct MultiChromFmIndex {
    /// Per-chromosome FM-Index structures.
    pub indices: Vec<FmIndex>,
}

impl MultiChromFmIndex {
    /// Creates an empty [`MultiChromFmIndex`].
    pub fn new() -> Self {
        Self {
            indices: Vec::new(),
        }
    }

    /// Adds a chromosome sequence to the assembly FM-Index.
    pub fn add_chromosome(&mut self, chr: &str, sequence: &[u8]) {
        self.indices.push(FmIndex::build(chr, sequence));
    }

    /// Locates all approximate matches `(chr, position, mismatches)` across all indexed chromosomes.
    pub fn locate_approximate_all(
        &self,
        pattern: &[u8],
        max_mismatches: usize,
    ) -> Vec<(String, u64, usize)> {
        let mut results = Vec::new();
        for idx in &self.indices {
            for (pos, mm) in idx.locate_approximate(pattern, max_mismatches) {
                results.push((idx.chr.clone(), pos, mm));
            }
        }
        results.sort_by(|a, b| a.2.cmp(&b.2).then_with(|| a.0.cmp(&b.0)).then_with(|| a.1.cmp(&b.1)));
        results
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fm_index_exact_and_bwt_properties() {
        let seq = b"GATGCGATCGATCGATCGATCGGGCATCGATCGATC";
        let fm = FmIndex::build("chr1", seq);
        assert_eq!(fm.bwt_ascii().len(), seq.len() + 1);

        // Exact matches of "GATCG" occur at positions 5, 9, 13, 17, 28
        let positions = fm.locate_exact(b"GATCG");
        assert_eq!(positions, vec![5, 9, 13, 17, 28]);
        assert_eq!(fm.count_exact(b"GATCG"), 5);
        assert_eq!(fm.count_exact(b"GGGCAT"), 1);
        assert_eq!(fm.locate_exact(b"GGGCAT"), vec![21]);
        assert_eq!(fm.count_exact(b"TTTTTT"), 0);
    }

    #[test]
    fn test_fm_index_approximate_backtracking_matches_linear_scan() {
        let seq = b"CCGTAGCTAGCTAGCTAGCTCGGGTAGCTAGCTAGAAATTTCCCGGG";
        let query = b"CCGTAGCTAGCTAGCTAGCT";
        let fm = FmIndex::build("chr9", seq);

        let hits = fm.locate_approximate(query, 2);
        // Exact match at 0 must be found with 0 mismatches
        assert!(hits.contains(&(0, 0)));

        // Verify every hit returned by FM-Index backtracking has <= 2 mismatches against seq
        for (pos, mm) in &hits {
            let p = *pos as usize;
            let actual_mm = seq[p..p + query.len()]
                .iter()
                .zip(query.iter())
                .filter(|(a, b)| a != b)
                .count();
            assert_eq!(actual_mm, *mm);
            assert!(*mm <= 2);
        }
    }
}
