//! Taxonomic classification types.

use serde::{Deserialize, Serialize};

/// Full eight-rank lineage of an organism.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Lineage {
    /// Cellular domain (e.g. "Eukaryota").
    pub domain: String,
    /// Kingdom (e.g. "Animalia").
    pub kingdom: String,
    /// Phylum (e.g. "Chordata").
    pub phylum: String,
    /// Class (e.g. "Mammalia").
    pub class: String,
    /// Order (e.g. "Primates").
    pub order: String,
    /// Family (e.g. "Hominidae").
    pub family: String,
    /// Genus (e.g. "Homo").
    pub genus: String,
    /// Species epithet (e.g. "sapiens").
    pub species: String,
}

/// One of the eight canonical taxonomic ranks.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TaxRank {
    /// Cellular domain.
    Domain,
    /// Kingdom rank.
    Kingdom,
    /// Phylum rank.
    Phylum,
    /// Class rank.
    Class,
    /// Order rank.
    Order,
    /// Family rank.
    Family,
    /// Genus rank.
    Genus,
    /// Species rank.
    Species,
}

/// A biological organism record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Organism {
    /// NCBI Taxonomy ID (e.g. 9606 for *Homo sapiens*).
    pub tax_id: u32,
    /// Latin binomial name (e.g. "Homo sapiens").
    pub scientific_name: String,
    /// Optional vernacular name (e.g. "Human").
    pub common_name: Option<String>,
    /// Eight-rank lineage.
    pub lineage: Lineage,
}

impl Organism {
    /// Returns `true` when this organism is *Homo sapiens* (tax_id 9606).
    pub fn is_human(&self) -> bool {
        self.tax_id == 9606
    }
}
