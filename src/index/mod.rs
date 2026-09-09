//! Genomic index structures for fast coordinate and sequence lookup.

pub mod interval_idx;
pub mod kmer_idx;

pub use interval_idx::{IndexEntry, IntervalIndex};
pub use kmer_idx::KmerIndex;
