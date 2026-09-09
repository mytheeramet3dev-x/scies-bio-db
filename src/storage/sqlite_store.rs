//! SQLite-backed persistent store for structured biological metadata.

use std::path::Path;

use rusqlite::{params, Connection, OptionalExtension};

use crate::model::gene::{Biotype, Gene, Transcript, TranscriptBiotype};
use crate::model::genome::{AssemblyStatus, Chromosome, ReferenceGenome, Strand};
use crate::model::taxonomy::{Lineage, Organism};
use crate::Result;

// ---------------------------------------------------------------------------
// CREATE TABLE strings (used in `init_schema`)
// ---------------------------------------------------------------------------

const SQL_CREATE_ORGANISMS: &str = "
CREATE TABLE IF NOT EXISTS organisms (
    tax_id          INTEGER PRIMARY KEY,
    scientific_name TEXT    NOT NULL,
    common_name     TEXT,
    domain          TEXT    NOT NULL,
    kingdom         TEXT    NOT NULL,
    phylum          TEXT    NOT NULL,
    class           TEXT    NOT NULL,
    ord             TEXT    NOT NULL,
    family          TEXT    NOT NULL,
    genus           TEXT    NOT NULL,
    species         TEXT    NOT NULL
);";

const SQL_CREATE_REFERENCE_GENOMES: &str = "
CREATE TABLE IF NOT EXISTS reference_genomes (
    assembly_name TEXT    PRIMARY KEY,
    tax_id        INTEGER NOT NULL REFERENCES organisms(tax_id),
    status        TEXT    NOT NULL,
    release_year  INTEGER NOT NULL
);";

const SQL_CREATE_CHROMOSOMES: &str = "
CREATE TABLE IF NOT EXISTS chromosomes (
    name              TEXT    NOT NULL,
    assembly_name     TEXT    NOT NULL REFERENCES reference_genomes(assembly_name),
    length            INTEGER NOT NULL,
    is_mitochondrial  INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (name, assembly_name)
);";

const SQL_CREATE_GENES: &str = "
CREATE TABLE IF NOT EXISTS genes (
    gene_id  TEXT    PRIMARY KEY,
    symbol   TEXT    NOT NULL,
    name     TEXT    NOT NULL,
    biotype  TEXT    NOT NULL,
    chr      TEXT    NOT NULL,
    start    INTEGER NOT NULL,
    end      INTEGER NOT NULL,
    strand   TEXT    NOT NULL,
    tax_id   INTEGER NOT NULL REFERENCES organisms(tax_id)
);";

const SQL_CREATE_TRANSCRIPTS: &str = "
CREATE TABLE IF NOT EXISTS transcripts (
    transcript_id TEXT    PRIMARY KEY,
    gene_id       TEXT    NOT NULL REFERENCES genes(gene_id),
    biotype       TEXT    NOT NULL,
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        TEXT    NOT NULL,
    exon_count    INTEGER NOT NULL
);";

const SQL_CREATE_EXONS: &str = "
CREATE TABLE IF NOT EXISTS exons (
    exon_id       TEXT    NOT NULL,
    transcript_id TEXT    NOT NULL REFERENCES transcripts(transcript_id),
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        TEXT    NOT NULL,
    exon_number   INTEGER NOT NULL,
    PRIMARY KEY (exon_id, transcript_id)
);";

const SQL_CREATE_PROTEIN_ENTRIES: &str = "
CREATE TABLE IF NOT EXISTS protein_entries (
    uniprot_accession TEXT    PRIMARY KEY,
    entry_name        TEXT    NOT NULL,
    gene_symbol       TEXT    NOT NULL,
    tax_id            INTEGER NOT NULL REFERENCES organisms(tax_id),
    status            TEXT    NOT NULL,
    sequence          TEXT    NOT NULL,
    length            INTEGER NOT NULL
);";

const SQL_CREATE_PROTEIN_DOMAINS: &str = "
CREATE TABLE IF NOT EXISTS protein_domains (
    domain_id         TEXT    NOT NULL,
    name              TEXT    NOT NULL,
    db                TEXT    NOT NULL,
    uniprot_accession TEXT    NOT NULL REFERENCES protein_entries(uniprot_accession),
    start_aa          INTEGER NOT NULL,
    end_aa            INTEGER NOT NULL,
    PRIMARY KEY (domain_id, uniprot_accession)
);";

const SQL_CREATE_VARIANT_RECORDS: &str = "
CREATE TABLE IF NOT EXISTS variant_records (
    variant_id   TEXT    PRIMARY KEY,
    chr          TEXT    NOT NULL,
    position     INTEGER NOT NULL,
    reference    TEXT    NOT NULL,
    alternate    TEXT    NOT NULL,
    variant_type TEXT    NOT NULL,
    tax_id       INTEGER NOT NULL REFERENCES organisms(tax_id)
);";

const SQL_CREATE_CLINICAL_VARIANTS: &str = "
CREATE TABLE IF NOT EXISTS clinical_variants (
    variant_id              TEXT    NOT NULL REFERENCES variant_records(variant_id),
    clinical_significance   TEXT    NOT NULL,
    condition               TEXT    NOT NULL,
    review_status           TEXT    NOT NULL,
    gene_symbol             TEXT,
    PRIMARY KEY (variant_id, condition)
);";

// ---------------------------------------------------------------------------
// SqliteStore
// ---------------------------------------------------------------------------

/// A thin wrapper around a [`rusqlite::Connection`] providing schema
/// initialisation and transactional access.
pub struct SqliteStore {
    /// The underlying SQLite connection.
    pub(crate) conn: std::sync::Mutex<Connection>,
}

impl SqliteStore {
    /// Opens (or creates) the SQLite database at `path`.
    pub fn open(path: &Path) -> Result<Self> {
        let conn = Connection::open(path)?;
        // Enable WAL for better concurrent read performance.
        conn.execute_batch("PRAGMA journal_mode=WAL; PRAGMA foreign_keys=ON;")?;
        Ok(Self {
            conn: std::sync::Mutex::new(conn),
        })
    }

    /// Creates all required tables if they do not already exist.
    ///
    /// This is idempotent and safe to call on every start-up.
    pub fn init_schema(&self) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let schema_batch = format!(
            "{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}",
            SQL_CREATE_ORGANISMS,
            SQL_CREATE_REFERENCE_GENOMES,
            SQL_CREATE_CHROMOSOMES,
            SQL_CREATE_GENES,
            SQL_CREATE_TRANSCRIPTS,
            SQL_CREATE_EXONS,
            SQL_CREATE_PROTEIN_ENTRIES,
            SQL_CREATE_PROTEIN_DOMAINS,
            SQL_CREATE_VARIANT_RECORDS,
            SQL_CREATE_CLINICAL_VARIANTS,
        );
        conn.execute_batch(&schema_batch)?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Organisms & Taxonomy
    // -----------------------------------------------------------------------

    pub fn insert_organism(&self, org: &Organism) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT OR REPLACE INTO organisms (tax_id, scientific_name, common_name, domain, kingdom, phylum, class, ord, family, genus, species)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            params![
                org.tax_id,
                org.scientific_name,
                org.common_name,
                org.lineage.domain,
                org.lineage.kingdom,
                org.lineage.phylum,
                org.lineage.class,
                org.lineage.order,
                org.lineage.family,
                org.lineage.genus,
                org.lineage.species,
            ],
        )?;
        Ok(())
    }

    pub fn get_organism(&self, tax_id: u32) -> Result<Option<Organism>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, scientific_name, common_name, domain, kingdom, phylum, class, ord, family, genus, species FROM organisms WHERE tax_id = ?1",
        )?;
        let org = stmt
            .query_row(params![tax_id], |row| {
                Ok(Organism {
                    tax_id: row.get(0)?,
                    scientific_name: row.get(1)?,
                    common_name: row.get(2)?,
                    lineage: Lineage {
                        domain: row.get(3)?,
                        kingdom: row.get(4)?,
                        phylum: row.get(5)?,
                        class: row.get(6)?,
                        order: row.get(7)?,
                        family: row.get(8)?,
                        genus: row.get(9)?,
                        species: row.get(10)?,
                    },
                })
            })
            .optional()?;
        Ok(org)
    }

    pub fn list_organisms(&self) -> Result<Vec<Organism>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, scientific_name, common_name, domain, kingdom, phylum, class, ord, family, genus, species FROM organisms ORDER BY tax_id",
        )?;
        let rows = stmt.query_map([], |row| {
            Ok(Organism {
                tax_id: row.get(0)?,
                scientific_name: row.get(1)?,
                common_name: row.get(2)?,
                lineage: Lineage {
                    domain: row.get(3)?,
                    kingdom: row.get(4)?,
                    phylum: row.get(5)?,
                    class: row.get(6)?,
                    order: row.get(7)?,
                    family: row.get(8)?,
                    genus: row.get(9)?,
                    species: row.get(10)?,
                },
            })
        })?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    // -----------------------------------------------------------------------
    // Reference Genomes & Chromosomes
    // -----------------------------------------------------------------------

    pub fn insert_reference_genome(&self, genome: &ReferenceGenome) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let status_str = match genome.status {
            AssemblyStatus::Primary => "Primary",
            AssemblyStatus::Alternate => "Alternate",
            AssemblyStatus::Patch => "Patch",
        };
        conn.execute(
            "INSERT OR REPLACE INTO reference_genomes (assembly_name, tax_id, status, release_year)
             VALUES (?1, ?2, ?3, ?4)",
            params![
                genome.assembly_name,
                genome.tax_id,
                status_str,
                genome.release_year
            ],
        )?;
        Ok(())
    }

    pub fn insert_chromosome(&self, chrom: &Chromosome) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT OR REPLACE INTO chromosomes (name, assembly_name, length, is_mitochondrial)
             VALUES (?1, ?2, ?3, ?4)",
            params![
                chrom.name,
                chrom.assembly_name,
                chrom.length,
                if chrom.is_mitochondrial { 1 } else { 0 }
            ],
        )?;
        Ok(())
    }

    pub fn list_chromosomes(&self, assembly: &str) -> Result<Vec<Chromosome>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT name, assembly_name, length, is_mitochondrial FROM chromosomes WHERE assembly_name = ?1 ORDER BY name",
        )?;
        let rows = stmt.query_map(params![assembly], |row| {
            let is_mito_int: i32 = row.get(3)?;
            Ok(Chromosome {
                name: row.get(0)?,
                assembly_name: row.get(1)?,
                length: row.get(2)?,
                is_mitochondrial: is_mito_int != 0,
            })
        })?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    // -----------------------------------------------------------------------
    // Genes & Transcripts
    // -----------------------------------------------------------------------

    pub fn insert_gene(&self, gene: &Gene) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let biotype_str = match &gene.biotype {
            Biotype::ProteinCoding => "protein_coding".to_string(),
            Biotype::LongNonCodingRna => "lncRNA".to_string(),
            Biotype::MicroRna => "miRNA".to_string(),
            Biotype::SmallNuclearRna => "snRNA".to_string(),
            Biotype::Pseudogene => "pseudogene".to_string(),
            Biotype::RibosomalRna => "rRNA".to_string(),
            Biotype::Other(s) => s.clone(),
        };
        let strand_str = match gene.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
        };
        conn.execute(
            "INSERT OR REPLACE INTO genes (gene_id, symbol, name, biotype, chr, start, end, strand, tax_id)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
            params![
                gene.gene_id,
                gene.symbol,
                gene.name,
                biotype_str,
                gene.chr,
                gene.start,
                gene.end,
                strand_str,
                gene.tax_id
            ],
        )?;
        Ok(())
    }

    pub fn find_gene_by_symbol(&self, symbol: &str) -> Result<Option<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, symbol, name, biotype, chr, start, end, strand, tax_id FROM genes WHERE symbol = ?1",
        )?;
        let gene = stmt
            .query_row(params![symbol], Self::row_to_gene)
            .optional()?;
        Ok(gene)
    }

    pub fn find_gene_by_id(&self, gene_id: &str) -> Result<Option<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, symbol, name, biotype, chr, start, end, strand, tax_id FROM genes WHERE gene_id = ?1",
        )?;
        let gene = stmt
            .query_row(params![gene_id], Self::row_to_gene)
            .optional()?;
        Ok(gene)
    }

    pub fn find_genes_by_region(&self, chr: &str, start: u64, end: u64) -> Result<Vec<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, symbol, name, biotype, chr, start, end, strand, tax_id FROM genes WHERE chr = ?1 AND start < ?3 AND end > ?2 ORDER BY start",
        )?;
        let rows = stmt.query_map(params![chr, start, end], Self::row_to_gene)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    fn row_to_gene(row: &rusqlite::Row) -> rusqlite::Result<Gene> {
        let biotype_raw: String = row.get(3)?;
        let biotype = match biotype_raw.as_str() {
            "protein_coding" => Biotype::ProteinCoding,
            "lncRNA" => Biotype::LongNonCodingRna,
            "miRNA" => Biotype::MicroRna,
            "snRNA" => Biotype::SmallNuclearRna,
            "pseudogene" => Biotype::Pseudogene,
            "rRNA" => Biotype::RibosomalRna,
            other => Biotype::Other(other.to_string()),
        };
        let strand_raw: String = row.get(7)?;
        let strand = if strand_raw == "-" {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        Ok(Gene {
            gene_id: row.get(0)?,
            symbol: row.get(1)?,
            name: row.get(2)?,
            biotype,
            chr: row.get(4)?,
            start: row.get(5)?,
            end: row.get(6)?,
            strand,
            tax_id: row.get(8)?,
        })
    }

    pub fn insert_transcript(&self, transcript: &Transcript) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let biotype_str = match &transcript.biotype {
            TranscriptBiotype::MRNA => "mRNA".to_string(),
            TranscriptBiotype::NCRNA => "ncRNA".to_string(),
            TranscriptBiotype::Pseudogene => "pseudogene".to_string(),
            TranscriptBiotype::Other(s) => s.clone(),
        };
        let strand_str = match transcript.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
        };
        conn.execute(
            "INSERT OR REPLACE INTO transcripts (transcript_id, gene_id, biotype, chr, start, end, strand, exon_count)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            params![
                transcript.transcript_id,
                transcript.gene_id,
                biotype_str,
                transcript.chr,
                transcript.start,
                transcript.end,
                strand_str,
                transcript.exon_count
            ],
        )?;
        Ok(())
    }

    pub fn find_transcripts_by_gene(&self, gene_id: &str) -> Result<Vec<Transcript>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT transcript_id, gene_id, biotype, chr, start, end, strand, exon_count FROM transcripts WHERE gene_id = ?1 ORDER BY start",
        )?;
        let rows = stmt.query_map(params![gene_id], |row| {
            let biotype_raw: String = row.get(2)?;
            let biotype = match biotype_raw.as_str() {
                "mRNA" => TranscriptBiotype::MRNA,
                "ncRNA" => TranscriptBiotype::NCRNA,
                "pseudogene" => TranscriptBiotype::Pseudogene,
                other => TranscriptBiotype::Other(other.to_string()),
            };
            let strand_raw: String = row.get(6)?;
            let strand = if strand_raw == "-" {
                Strand::Reverse
            } else {
                Strand::Forward
            };
            Ok(Transcript {
                transcript_id: row.get(0)?,
                gene_id: row.get(1)?,
                biotype,
                chr: row.get(3)?,
                start: row.get(4)?,
                end: row.get(5)?,
                strand,
                exon_count: row.get(7)?,
            })
        })?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_sqlite_init_and_crud() {
        let dir = tempdir().unwrap();
        let db_path = dir.path().join("test.db");
        let store = SqliteStore::open(&db_path).unwrap();
        store.init_schema().unwrap();

        let org = Organism {
            tax_id: 9606,
            scientific_name: "Homo sapiens".to_string(),
            common_name: Some("Human".to_string()),
            lineage: Lineage {
                domain: "Eukaryota".to_string(),
                kingdom: "Animalia".to_string(),
                phylum: "Chordata".to_string(),
                class: "Mammalia".to_string(),
                order: "Primates".to_string(),
                family: "Hominidae".to_string(),
                genus: "Homo".to_string(),
                species: "sapiens".to_string(),
            },
        };
        store.insert_organism(&org).unwrap();

        let retrieved = store.get_organism(9606).unwrap().unwrap();
        assert_eq!(retrieved.scientific_name, "Homo sapiens");
        assert_eq!(retrieved.common_name, Some("Human".to_string()));
        assert_eq!(retrieved.lineage.genus, "Homo");
    }
}
