//! SQLite-backed persistent store for structured biological metadata.

use std::path::Path;

use rusqlite::{params, Connection, OptionalExtension};

use crate::model::gene::{Biotype, Exon, Gene, Transcript, TranscriptBiotype};
use crate::model::genome::{AssemblyStatus, Chromosome, ReferenceGenome, Strand};
use crate::model::protein::{Domain, ProteinEntry, ReviewStatus};
use crate::model::taxonomy::{Lineage, Organism};
use crate::model::variant::{ClinicalSignificance, ClinicalVariant, VariantRecord, VariantType};
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
    gene_id       TEXT    NOT NULL,
    assembly_name TEXT    NOT NULL REFERENCES reference_genomes(assembly_name),
    symbol        TEXT    NOT NULL,
    name          TEXT    NOT NULL,
    biotype       TEXT    NOT NULL,
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        TEXT    NOT NULL,
    tax_id        INTEGER NOT NULL REFERENCES organisms(tax_id),
    PRIMARY KEY (gene_id, assembly_name)
);
CREATE INDEX IF NOT EXISTS idx_genes_region ON genes(tax_id, assembly_name, chr, start, end);
CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(tax_id, symbol);
CREATE INDEX IF NOT EXISTS idx_genes_symbol_only ON genes(symbol);
";

const SQL_CREATE_TRANSCRIPTS: &str = "
CREATE TABLE IF NOT EXISTS transcripts (
    transcript_id TEXT    NOT NULL,
    gene_id       TEXT    NOT NULL,
    assembly_name TEXT    NOT NULL,
    biotype       TEXT    NOT NULL,
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        TEXT    NOT NULL,
    exon_count    INTEGER NOT NULL,
    PRIMARY KEY (transcript_id, assembly_name),
    FOREIGN KEY (gene_id, assembly_name) REFERENCES genes(gene_id, assembly_name) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_transcripts_gene ON transcripts(gene_id, assembly_name);
";

const SQL_CREATE_EXONS: &str = "
CREATE TABLE IF NOT EXISTS exons (
    exon_id       TEXT    NOT NULL,
    transcript_id TEXT    NOT NULL,
    assembly_name TEXT    NOT NULL,
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        TEXT    NOT NULL,
    exon_number   INTEGER NOT NULL,
    PRIMARY KEY (exon_id, transcript_id, assembly_name),
    FOREIGN KEY (transcript_id, assembly_name) REFERENCES transcripts(transcript_id, assembly_name) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_exons_transcript ON exons(transcript_id, assembly_name);
";

const SQL_CREATE_PROTEIN_ENTRIES: &str = "
CREATE TABLE IF NOT EXISTS protein_entries (
    uniprot_accession TEXT    PRIMARY KEY,
    entry_name        TEXT    NOT NULL,
    gene_symbol       TEXT    NOT NULL,
    tax_id            INTEGER NOT NULL REFERENCES organisms(tax_id),
    status            TEXT    NOT NULL,
    sequence          TEXT    NOT NULL,
    length            INTEGER NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_proteins_gene ON protein_entries(gene_symbol);
";

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
    variant_id    TEXT    NOT NULL,
    assembly_name TEXT    NOT NULL REFERENCES reference_genomes(assembly_name),
    chr           TEXT    NOT NULL,
    position      INTEGER NOT NULL,
    reference     TEXT    NOT NULL,
    alternate     TEXT    NOT NULL,
    variant_type  TEXT    NOT NULL,
    tax_id        INTEGER NOT NULL REFERENCES organisms(tax_id),
    PRIMARY KEY (variant_id, assembly_name)
);
CREATE INDEX IF NOT EXISTS idx_variants_region ON variant_records(tax_id, assembly_name, chr, position);
";

const SQL_CREATE_CLINICAL_VARIANTS: &str = "
CREATE TABLE IF NOT EXISTS clinical_variants (
    variant_id              TEXT    NOT NULL,
    condition               TEXT    NOT NULL,
    clinical_significance   TEXT    NOT NULL,
    review_status           TEXT    NOT NULL,
    gene_symbol             TEXT,
    PRIMARY KEY (variant_id, condition)
);
CREATE INDEX IF NOT EXISTS idx_clinical_gene ON clinical_variants(gene_symbol);
";

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
        // Enable WAL, foreign keys, and 5-second busy timeout.
        conn.execute_batch(
            "PRAGMA journal_mode=WAL; PRAGMA foreign_keys=ON; PRAGMA busy_timeout=5000;",
        )?;
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
            "INSERT INTO organisms (tax_id, scientific_name, common_name, domain, kingdom, phylum, class, ord, family, genus, species)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)
             ON CONFLICT (tax_id) DO UPDATE SET
                scientific_name = excluded.scientific_name,
                common_name = excluded.common_name,
                domain = excluded.domain,
                kingdom = excluded.kingdom,
                phylum = excluded.phylum,
                class = excluded.class,
                ord = excluded.ord,
                family = excluded.family,
                genus = excluded.genus,
                species = excluded.species",
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
            "INSERT INTO reference_genomes (assembly_name, tax_id, status, release_year)
             VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT (assembly_name) DO UPDATE SET
                tax_id = excluded.tax_id,
                status = excluded.status,
                release_year = excluded.release_year",
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
            "INSERT INTO chromosomes (name, assembly_name, length, is_mitochondrial)
             VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT (name, assembly_name) DO UPDATE SET
                length = excluded.length,
                is_mitochondrial = excluded.is_mitochondrial",
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
    // Genes, Transcripts & Exons
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
            "INSERT INTO genes (gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)
             ON CONFLICT (gene_id, assembly_name) DO UPDATE SET
                symbol = excluded.symbol,
                name = excluded.name,
                biotype = excluded.biotype,
                chr = excluded.chr,
                start = excluded.start,
                end = excluded.end,
                strand = excluded.strand,
                tax_id = excluded.tax_id",
            params![
                gene.gene_id,
                gene.assembly_name,
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

    /// Finds a gene by tax_id and HGNC symbol (species-isolated).
    pub fn find_gene_by_symbol(&self, tax_id: u32, symbol: &str) -> Result<Option<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id FROM genes WHERE tax_id = ?1 AND symbol = ?2 LIMIT 1",
        )?;
        let gene = stmt
            .query_row(params![tax_id, symbol], Self::row_to_gene)
            .optional()?;
        Ok(gene)
    }

    /// Finds a gene by symbol across all species (returns first match).
    pub fn find_gene_by_symbol_any(&self, symbol: &str) -> Result<Option<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id FROM genes WHERE symbol = ?1 LIMIT 1",
        )?;
        let gene = stmt
            .query_row(params![symbol], Self::row_to_gene)
            .optional()?;
        Ok(gene)
    }

    /// Finds a gene by identifier across all assemblies (returns first match).
    pub fn find_gene_by_id(&self, gene_id: &str) -> Result<Option<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id FROM genes WHERE gene_id = ?1 LIMIT 1",
        )?;
        let gene = stmt
            .query_row(params![gene_id], Self::row_to_gene)
            .optional()?;
        Ok(gene)
    }

    /// Finds a gene by identifier scoped to a specific assembly.
    pub fn find_gene_by_id_scoped(&self, gene_id: &str, assembly: &str) -> Result<Option<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id FROM genes WHERE gene_id = ?1 AND assembly_name = ?2",
        )?;
        let gene = stmt
            .query_row(params![gene_id, assembly], Self::row_to_gene)
            .optional()?;
        Ok(gene)
    }

    /// Finds genes overlapping interval [start, end) on chromosome `chr`, strictly scoped by `tax_id` and `assembly`.
    pub fn find_genes_by_region(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id
             FROM genes
             WHERE tax_id = ?1 AND assembly_name = ?2 AND chr = ?3 AND start < ?5 AND end > ?4
             ORDER BY start",
        )?;
        let rows = stmt.query_map(
            params![tax_id, assembly, chr, start, end],
            Self::row_to_gene,
        )?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    /// Unscoped region query (convenience fallback).
    pub fn find_genes_by_region_unscoped(
        &self,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id
             FROM genes
             WHERE chr = ?1 AND start < ?3 AND end > ?2
             ORDER BY start",
        )?;
        let rows = stmt.query_map(params![chr, start, end], Self::row_to_gene)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    fn row_to_gene(row: &rusqlite::Row) -> rusqlite::Result<Gene> {
        let biotype_raw: String = row.get(4)?;
        let biotype = match biotype_raw.as_str() {
            "protein_coding" => Biotype::ProteinCoding,
            "lncRNA" => Biotype::LongNonCodingRna,
            "miRNA" => Biotype::MicroRna,
            "snRNA" => Biotype::SmallNuclearRna,
            "pseudogene" => Biotype::Pseudogene,
            "rRNA" => Biotype::RibosomalRna,
            other => Biotype::Other(other.to_string()),
        };
        let strand_raw: String = row.get(8)?;
        let strand = if strand_raw == "-" {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        Ok(Gene {
            gene_id: row.get(0)?,
            assembly_name: row.get(1)?,
            symbol: row.get(2)?,
            name: row.get(3)?,
            biotype,
            chr: row.get(5)?,
            start: row.get(6)?,
            end: row.get(7)?,
            strand,
            tax_id: row.get(9)?,
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
            "INSERT INTO transcripts (transcript_id, gene_id, assembly_name, biotype, chr, start, end, strand, exon_count)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)
             ON CONFLICT (transcript_id, assembly_name) DO UPDATE SET
                gene_id = excluded.gene_id,
                biotype = excluded.biotype,
                chr = excluded.chr,
                start = excluded.start,
                end = excluded.end,
                strand = excluded.strand,
                exon_count = excluded.exon_count",
            params![
                transcript.transcript_id,
                transcript.gene_id,
                transcript.assembly_name,
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

    pub fn insert_exon(&self, exon: &Exon) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let strand_str = match exon.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
        };
        conn.execute(
            "INSERT INTO exons (exon_id, transcript_id, assembly_name, chr, start, end, strand, exon_number)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)
             ON CONFLICT (exon_id, transcript_id, assembly_name) DO UPDATE SET
                chr = excluded.chr,
                start = excluded.start,
                end = excluded.end,
                strand = excluded.strand,
                exon_number = excluded.exon_number",
            params![
                exon.exon_id,
                exon.transcript_id,
                exon.assembly_name,
                exon.chr,
                exon.start,
                exon.end,
                strand_str,
                exon.exon_number,
            ],
        )?;
        Ok(())
    }

    pub fn find_transcripts_by_gene(&self, gene_id: &str) -> Result<Vec<Transcript>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT transcript_id, gene_id, assembly_name, biotype, chr, start, end, strand, exon_count
             FROM transcripts WHERE gene_id = ?1 ORDER BY start",
        )?;
        let rows = stmt.query_map(params![gene_id], Self::row_to_transcript)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    pub fn find_transcripts_by_gene_scoped(
        &self,
        gene_id: &str,
        assembly: &str,
    ) -> Result<Vec<Transcript>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT transcript_id, gene_id, assembly_name, biotype, chr, start, end, strand, exon_count
             FROM transcripts WHERE gene_id = ?1 AND assembly_name = ?2 ORDER BY start",
        )?;
        let rows = stmt.query_map(params![gene_id, assembly], Self::row_to_transcript)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    fn row_to_transcript(row: &rusqlite::Row) -> rusqlite::Result<Transcript> {
        let biotype_raw: String = row.get(3)?;
        let biotype = match biotype_raw.as_str() {
            "mRNA" => TranscriptBiotype::MRNA,
            "ncRNA" => TranscriptBiotype::NCRNA,
            "pseudogene" => TranscriptBiotype::Pseudogene,
            other => TranscriptBiotype::Other(other.to_string()),
        };
        let strand_raw: String = row.get(7)?;
        let strand = if strand_raw == "-" {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        Ok(Transcript {
            transcript_id: row.get(0)?,
            gene_id: row.get(1)?,
            assembly_name: row.get(2)?,
            biotype,
            chr: row.get(4)?,
            start: row.get(5)?,
            end: row.get(6)?,
            strand,
            exon_count: row.get(8)?,
        })
    }

    // -----------------------------------------------------------------------
    // Protein Entries & Domains
    // -----------------------------------------------------------------------

    pub fn insert_protein_entry(&self, entry: &ProteinEntry) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let status_str = match entry.status {
            ReviewStatus::SwissProt => "SwissProt",
            ReviewStatus::TrEMBL => "TrEMBL",
        };
        conn.execute(
            "INSERT INTO protein_entries (uniprot_accession, entry_name, gene_symbol, tax_id, status, sequence, length)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)
             ON CONFLICT (uniprot_accession) DO UPDATE SET
                entry_name = excluded.entry_name,
                gene_symbol = excluded.gene_symbol,
                tax_id = excluded.tax_id,
                status = excluded.status,
                sequence = excluded.sequence,
                length = excluded.length",
            params![
                entry.uniprot_accession,
                entry.entry_name,
                entry.gene_symbol,
                entry.tax_id,
                status_str,
                entry.sequence,
                entry.length,
            ],
        )?;
        Ok(())
    }

    pub fn insert_protein_domain(&self, domain: &Domain) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO protein_domains (domain_id, name, db, uniprot_accession, start_aa, end_aa)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)
             ON CONFLICT (domain_id, uniprot_accession) DO UPDATE SET
                name = excluded.name,
                db = excluded.db,
                start_aa = excluded.start_aa,
                end_aa = excluded.end_aa",
            params![
                domain.domain_id,
                domain.name,
                domain.db,
                domain.uniprot_accession,
                domain.start_aa,
                domain.end_aa,
            ],
        )?;
        Ok(())
    }

    pub fn find_protein_by_accession(&self, accession: &str) -> Result<Option<ProteinEntry>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT uniprot_accession, entry_name, gene_symbol, tax_id, status, sequence, length
             FROM protein_entries WHERE uniprot_accession = ?1",
        )?;
        let entry = stmt
            .query_row(params![accession], Self::row_to_protein)
            .optional()?;
        Ok(entry)
    }

    pub fn find_proteins_by_gene(&self, gene_symbol: &str) -> Result<Vec<ProteinEntry>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT uniprot_accession, entry_name, gene_symbol, tax_id, status, sequence, length
             FROM protein_entries WHERE gene_symbol = ?1 ORDER BY uniprot_accession",
        )?;
        let rows = stmt.query_map(params![gene_symbol], Self::row_to_protein)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    pub fn find_domains_by_accession(&self, accession: &str) -> Result<Vec<Domain>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT domain_id, name, db, uniprot_accession, start_aa, end_aa
             FROM protein_domains WHERE uniprot_accession = ?1 ORDER BY start_aa",
        )?;
        let rows = stmt.query_map(params![accession], |row| {
            Ok(Domain {
                domain_id: row.get(0)?,
                name: row.get(1)?,
                db: row.get(2)?,
                uniprot_accession: row.get(3)?,
                start_aa: row.get(4)?,
                end_aa: row.get(5)?,
            })
        })?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    fn row_to_protein(row: &rusqlite::Row) -> rusqlite::Result<ProteinEntry> {
        let status_str: String = row.get(4)?;
        let status = match status_str.as_str() {
            "SwissProt" => ReviewStatus::SwissProt,
            _ => ReviewStatus::TrEMBL,
        };
        Ok(ProteinEntry {
            uniprot_accession: row.get(0)?,
            entry_name: row.get(1)?,
            gene_symbol: row.get(2)?,
            tax_id: row.get(3)?,
            status,
            sequence: row.get(5)?,
            length: row.get(6)?,
            subcellular_locations: Vec::new(),
        })
    }

    // -----------------------------------------------------------------------
    // Variants & Clinical Interpretations
    // -----------------------------------------------------------------------

    pub fn insert_variant_record(&self, variant: &VariantRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let type_str = match variant.variant_type {
            VariantType::Snp => "SNP",
            VariantType::Insertion => "Insertion",
            VariantType::Deletion => "Deletion",
            VariantType::Mnp => "MNP",
            VariantType::StructuralVariant => "SV",
        };
        conn.execute(
            "INSERT INTO variant_records (variant_id, assembly_name, chr, position, reference, alternate, variant_type, tax_id)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)
             ON CONFLICT (variant_id, assembly_name) DO UPDATE SET
                chr = excluded.chr,
                position = excluded.position,
                reference = excluded.reference,
                alternate = excluded.alternate,
                variant_type = excluded.variant_type,
                tax_id = excluded.tax_id",
            params![
                variant.variant_id,
                variant.assembly_name,
                variant.chr,
                variant.position,
                variant.reference,
                variant.alternate,
                type_str,
                variant.tax_id,
            ],
        )?;
        Ok(())
    }

    pub fn insert_clinical_variant(&self, clinical: &ClinicalVariant) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let sig_str = match clinical.clinical_significance {
            ClinicalSignificance::Pathogenic => "Pathogenic",
            ClinicalSignificance::LikelyPathogenic => "LikelyPathogenic",
            ClinicalSignificance::VUS => "VUS",
            ClinicalSignificance::LikelyBenign => "LikelyBenign",
            ClinicalSignificance::Benign => "Benign",
            ClinicalSignificance::Conflicting => "Conflicting",
            ClinicalSignificance::Unknown => "Unknown",
        };
        conn.execute(
            "INSERT INTO clinical_variants (variant_id, condition, clinical_significance, review_status, gene_symbol)
             VALUES (?1, ?2, ?3, ?4, ?5)
             ON CONFLICT (variant_id, condition) DO UPDATE SET
                clinical_significance = excluded.clinical_significance,
                review_status = excluded.review_status,
                gene_symbol = excluded.gene_symbol",
            params![
                clinical.variant_id,
                clinical.condition,
                sig_str,
                clinical.review_status,
                clinical.gene_symbol,
            ],
        )?;
        Ok(())
    }

    pub fn find_variant_by_id(&self, variant_id: &str) -> Result<Option<VariantRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT variant_id, assembly_name, chr, position, reference, alternate, variant_type, tax_id
             FROM variant_records WHERE variant_id = ?1 LIMIT 1",
        )?;
        let variant = stmt
            .query_row(params![variant_id], Self::row_to_variant)
            .optional()?;
        Ok(variant)
    }

    pub fn find_variant_by_id_scoped(
        &self,
        variant_id: &str,
        assembly: &str,
    ) -> Result<Option<VariantRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT variant_id, assembly_name, chr, position, reference, alternate, variant_type, tax_id
             FROM variant_records WHERE variant_id = ?1 AND assembly_name = ?2",
        )?;
        let variant = stmt
            .query_row(params![variant_id, assembly], Self::row_to_variant)
            .optional()?;
        Ok(variant)
    }

    pub fn find_variants_by_region(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<VariantRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT variant_id, assembly_name, chr, position, reference, alternate, variant_type, tax_id
             FROM variant_records
             WHERE tax_id = ?1 AND assembly_name = ?2 AND chr = ?3 AND position >= ?4 AND position < ?5
             ORDER BY position",
        )?;
        let rows = stmt.query_map(
            params![tax_id, assembly, chr, start, end],
            Self::row_to_variant,
        )?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    pub fn find_clinical_variants_by_gene(
        &self,
        gene_symbol: &str,
    ) -> Result<Vec<ClinicalVariant>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT variant_id, condition, clinical_significance, review_status, gene_symbol
             FROM clinical_variants WHERE gene_symbol = ?1 ORDER BY condition",
        )?;
        let rows = stmt.query_map(params![gene_symbol], Self::row_to_clinical_variant)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    fn row_to_variant(row: &rusqlite::Row) -> rusqlite::Result<VariantRecord> {
        let type_str: String = row.get(6)?;
        let variant_type = match type_str.as_str() {
            "Insertion" => VariantType::Insertion,
            "Deletion" => VariantType::Deletion,
            "MNP" => VariantType::Mnp,
            "SV" => VariantType::StructuralVariant,
            _ => VariantType::Snp,
        };
        Ok(VariantRecord {
            variant_id: row.get(0)?,
            assembly_name: row.get(1)?,
            chr: row.get(2)?,
            position: row.get(3)?,
            reference: row.get(4)?,
            alternate: row.get(5)?,
            variant_type,
            tax_id: row.get(7)?,
        })
    }

    fn row_to_clinical_variant(row: &rusqlite::Row) -> rusqlite::Result<ClinicalVariant> {
        let sig_str: String = row.get(2)?;
        let clinical_significance = match sig_str.as_str() {
            "Pathogenic" => ClinicalSignificance::Pathogenic,
            "LikelyPathogenic" => ClinicalSignificance::LikelyPathogenic,
            "VUS" => ClinicalSignificance::VUS,
            "LikelyBenign" => ClinicalSignificance::LikelyBenign,
            "Benign" => ClinicalSignificance::Benign,
            "Conflicting" => ClinicalSignificance::Conflicting,
            _ => ClinicalSignificance::Unknown,
        };
        Ok(ClinicalVariant {
            variant_id: row.get(0)?,
            condition: row.get(1)?,
            clinical_significance,
            review_status: row.get(3)?,
            gene_symbol: row.get(4)?,
        })
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
