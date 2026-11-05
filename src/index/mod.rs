//! Genomic index structures for fast coordinate and sequence lookup.

pub mod fm_idx;
pub mod interval_idx;
pub mod kmer_idx;

pub use fm_idx::{FmIndex, MultiChromFmIndex};
pub use interval_idx::{IndexEntry, IntervalIndex};
pub use kmer_idx::{KmerHit, KmerIndex};
