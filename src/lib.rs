//! # scies-bio-db
//!
//! A biological database engine for all living organisms, starting from *Homo sapiens*.
//!
//! Built on top of [`scies-bio-th`] for core data types and parsers.
//!
//! ## Quick Start
//!
//! ```no_run
//! use scies_bio_db::BioDb;
//! use std::path::Path;
//!
//! let db = BioDb::open(Path::new("./human_grch38")).unwrap();
//! let gene = db.genes().find_by_symbol("TP53").unwrap();
//! ```

pub mod api;
pub mod error;
pub mod index;
pub mod ingest;
pub mod model;
pub mod query;
pub mod storage;

pub use api::BioDb;
pub use error::{BioDbError, Result};
