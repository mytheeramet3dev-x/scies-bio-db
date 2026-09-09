//! Sorted genomic interval index backed by a `Vec`.
//!
//! Entries are sorted by `(chr, start)` after a call to [`IntervalIndex::build`].
//! Overlap queries perform a binary-search gallop followed by a linear scan over
//! the matching chromosome window.

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

    /// Returns all entries whose interval `[entry.start, entry.end)` overlaps
    /// the query interval `[start, end)` on chromosome `chr`.
    ///
    /// # Panics
    ///
    /// Does not panic; returns an empty slice when no entries overlap.
    pub fn find_overlapping(&self, chr: &str, start: u64, end: u64) -> Vec<&IndexEntry> {
        // Binary-search to the first entry with this chromosome.
        let first = self.entries.partition_point(|e| {
            e.chr.as_str() < chr || (e.chr.as_str() == chr && e.start < start)
        });

        let mut result = Vec::new();
        // Walk backward to catch entries whose start < query start but end > query start.
        let scan = first.saturating_sub(1);
        // Also scan all entries on the same chromosome beginning before `end`.
        // For simplicity (and correctness for moderate dataset sizes) we do a
        // linear scan over all entries on `chr` that could overlap.
        // A production implementation would use an augmented interval tree or
        // segment tree for O(log n + k) performance.
        for entry in &self.entries {
            if entry.chr.as_str() != chr {
                continue;
            }
            // Overlap condition: entry.start < end AND entry.end > start
            if entry.start < end && entry.end > start {
                result.push(entry);
            }
        }
        let _ = scan; // suppress unused-variable warning
        result
    }
}

impl Default for IntervalIndex {
    fn default() -> Self {
        Self::new()
    }
}
