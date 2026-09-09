//! Error types for `scies-bio-db`.

use thiserror::Error;

/// All errors that can arise within the biological database engine.
#[derive(Debug, Error)]
pub enum BioDbError {
    /// An underlying I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// A SQLite / rusqlite error.
    #[error("SQLite error: {0}")]
    Sql(#[from] rusqlite::Error),

    /// A biological-type error from `scies_bio_th`.
    #[error("Bio error: {0}")]
    Bio(#[from] scies_bio_th::BioError),

    /// A 2-bit encoding/decoding error.
    #[error("Encoding error: {0}")]
    Encoding(String),

    /// The requested entity was not found in the store.
    #[error("Not found: {0}")]
    NotFound(String),

    /// A genomic coordinate was out of range or otherwise invalid.
    #[error("Invalid coordinate: {chr}:{start}-{end}")]
    InvalidCoordinate { chr: String, start: u64, end: u64 },
}

/// Convenience `Result` alias used throughout the crate.
pub type Result<T> = std::result::Result<T, BioDbError>;
