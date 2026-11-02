//! Sorted genomic interval index backed by a `Vec`.
//!
//! Entries are sorted by `(chr, start)` after a call to [`IntervalIndex::build`].
//! Overlap queries perform a binary-search gallop followed by a linear scan over
//! the matching chromosome window.

use crate::model::gene::{Gene, Transcript};

/// A single indexed interval record associating a genomic coordinate range with
/// an entity identifier.
#[derive(Debug, Clone)]
pub struct IndexEntry {
    /// Chromosome / contig name.
    pub chr: String,
    /// 0-based start (inclusive).
    pub start: u64,
    /// 0-based end (exclusive).
    pub end: u64,
    /// Opaque identifier of the entity this interval belongs to (e.g. gene_id).
    pub payload_id: String,
}

/// A genomic interval index supporting overlap queries.
///
/// Entries must be inserted via [`IntervalIndex::insert`] and the index rebuilt
/// with [`IntervalIndex::build`] before queries are issued.
pub struct IntervalIndex {
    /// Sorted interval entries.
    entries: Vec<IndexEntry>,
}

impl IntervalIndex {
    /// Creates an empty [`IntervalIndex`].
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Appends `entry` to the index.
    ///
    /// The index must be rebuilt with [`Self::build`] after all inserts are
    /// complete before calling [`Self::find_overlapping`].
    pub fn insert(&mut self, entry: IndexEntry) {
        self.entries.push(entry);
    }

    /// Sorts all entries by `(chr, start)` to enable binary-search lookups.
    pub fn build(&mut self) {
        self.entries
            .sort_unstable_by(|a, b| a.chr.cmp(&b.chr).then(a.start.cmp(&b.start)));
    }

    /// Returns the number of entries in the index.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if the index contains no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Builds a sorted [`IntervalIndex`] from a slice of [`Gene`] records.
    pub fn from_genes(genes: &[Gene]) -> Self {
        let mut index = Self {
            entries: genes
                .iter()
                .map(|g| IndexEntry {
                    chr: g.chr.clone(),
                    start: g.start,
                    end: g.end,
                    payload_id: g.gene_id.clone(),
                })
                .collect(),
        };
        index.build();
        index
    }

    /// Builds a sorted [`IntervalIndex`] from a slice of [`Transcript`] records.
    pub fn from_transcripts(transcripts: &[Transcript]) -> Self {
        let mut index = Self {
            entries: transcripts
                .iter()
                .map(|t| IndexEntry {
                    chr: t.chr.clone(),
                    start: t.start,
                    end: t.end,
                    payload_id: t.transcript_id.clone(),
                })
                .collect(),
        };
        index.build();
        index
    }

    /// Returns all entries whose interval `[entry.start, entry.end)` overlaps
    /// the query interval `[start, end)` on chromosome `chr`.
    ///
    /// # Panics
    ///
    /// Does not panic; returns an empty vector when no entries overlap.
    pub fn find_overlapping(&self, chr: &str, start: u64, end: u64) -> Vec<&IndexEntry> {
        // Binary search to find the first entry belonging to `chr`
        let chr_start = self.entries.partition_point(|e| e.chr.as_str() < chr);
        let mut result = Vec::new();

        for entry in &self.entries[chr_start..] {
            if entry.chr.as_str() != chr {
                break; // Reached next chromosome in sorted list
            }
            if entry.start >= end {
                break; // Sorted by start: all subsequent entries start at or after query end
            }
            // Overlap condition: entry.start < end (checked above) AND entry.end > start
            if entry.end > start {
                result.push(entry);
            }
        }
        result
    }
}

impl Default for IntervalIndex {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_index_overlapping() {
        let mut index = IntervalIndex::new();
        index.insert(IndexEntry {
            chr: "chr1".to_string(),
            start: 100,
            end: 200,
            payload_id: "gene1".to_string(),
        });
        index.insert(IndexEntry {
            chr: "chr1".to_string(),
            start: 150,
            end: 250,
            payload_id: "gene2".to_string(),
        });
        index.insert(IndexEntry {
            chr: "chr1".to_string(),
            start: 300,
            end: 400,
            payload_id: "gene3".to_string(),
        });
        index.insert(IndexEntry {
            chr: "chr2".to_string(),
            start: 100,
            end: 200,
            payload_id: "gene4".to_string(),
        });

        index.build();

        // Overlaps gene1 and gene2
        let hits = index.find_overlapping("chr1", 120, 160);
        let ids: Vec<&str> = hits.iter().map(|e| e.payload_id.as_str()).collect();
        assert_eq!(ids, vec!["gene1", "gene2"]);

        // Isolated hit for gene3
        let hits3 = index.find_overlapping("chr1", 350, 450);
        assert_eq!(hits3.len(), 1);
        assert_eq!(hits3[0].payload_id, "gene3");

        // Chr2 isolation
        let hits_chr2 = index.find_overlapping("chr2", 50, 150);
        assert_eq!(hits_chr2.len(), 1);
        assert_eq!(hits_chr2[0].payload_id, "gene4");

        // No overlap
        let misses = index.find_overlapping("chr1", 500, 600);
        assert!(misses.is_empty());
    }
}
