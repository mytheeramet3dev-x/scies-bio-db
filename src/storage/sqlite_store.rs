//! SQLite-backed persistent store for structured biological metadata.

use std::path::Path;

use rusqlite::{params, Connection, OptionalExtension};

use crate::model::epigenome::{
    DatasetManifest, EpigeneticFeature, EpigeneticTrackType, ProteinStructureRecord,
};
use crate::model::gene::{Biotype, Exon, Gene, GeneFeature, Transcript, TranscriptBiotype};
use crate::model::genome::{Assembly, AssemblyStatus, Chromosome, ReferenceGenome, SequenceRef, Strand};
use crate::model::longevity::LongevityTarget;
use crate::model::protein::{Domain, ProteinEntry, ReviewStatus};
use crate::model::regulatory::{
    MultiOmicsTrackRecord, PharmaCompoundRecord, RegulatoryEdgeRecord, SingleCellAtlasRecord,
    SpatialCellRecord, SpatialTranscriptRecord, SyntheticGameteRecord,
};
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

const SQL_CREATE_ASSEMBLIES: &str = "
CREATE TABLE IF NOT EXISTS assemblies (
    tax_id        INTEGER NOT NULL,
    assembly_name TEXT    NOT NULL,
    source        TEXT,
    created_at    TEXT,
    notes         TEXT,
    PRIMARY KEY (tax_id, assembly_name)
);
CREATE INDEX IF NOT EXISTS idx_assemblies_tax ON assemblies(tax_id);
";

const SQL_CREATE_SEQUENCES: &str = "
CREATE TABLE IF NOT EXISTS sequences (
    tax_id        INTEGER NOT NULL,
    assembly_name TEXT    NOT NULL,
    seq_id        TEXT    NOT NULL,
    length_bp     INTEGER NOT NULL,
    path          TEXT    NOT NULL,
    sha256        TEXT    NOT NULL,
    source_uri    TEXT,
    PRIMARY KEY (tax_id, assembly_name, seq_id)
);
CREATE INDEX IF NOT EXISTS idx_sequences_lookup ON sequences(tax_id, assembly_name, seq_id);
";

const SQL_CREATE_GENES: &str = "
CREATE TABLE IF NOT EXISTS genes (
    gene_id       TEXT    NOT NULL,
    assembly_name TEXT    NOT NULL,
    symbol        TEXT    NOT NULL,
    name          TEXT    NOT NULL,
    biotype       TEXT    NOT NULL,
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        TEXT    NOT NULL,
    tax_id        INTEGER NOT NULL,
    feature_type  TEXT    NOT NULL DEFAULT 'gene',
    source        TEXT,
    PRIMARY KEY (gene_id, assembly_name)
);
CREATE INDEX IF NOT EXISTS idx_genes_region ON genes(tax_id, assembly_name, chr, start, end);
CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes(tax_id, assembly_name, symbol);
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

const SQL_CREATE_EPIGENETIC_FEATURES: &str = "
CREATE TABLE IF NOT EXISTS epigenetic_features (
    feature_id    TEXT    NOT NULL,
    tax_id        INTEGER NOT NULL REFERENCES organisms(tax_id),
    assembly_name TEXT    NOT NULL REFERENCES reference_genomes(assembly_name),
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    track_type    TEXT    NOT NULL,
    mark_label    TEXT    NOT NULL,
    signal_value  REAL    NOT NULL,
    PRIMARY KEY (feature_id, assembly_name)
);
CREATE INDEX IF NOT EXISTS idx_epi_region ON epigenetic_features(tax_id, assembly_name, chr, start, end);
";

const SQL_CREATE_PROTEIN_STRUCTURES: &str = "
CREATE TABLE IF NOT EXISTS protein_structures (
    structure_id        TEXT    PRIMARY KEY,
    gene_symbol         TEXT    NOT NULL,
    tax_id              INTEGER NOT NULL REFERENCES organisms(tax_id),
    resolution_angstrom REAL    NOT NULL,
    residue_count       INTEGER NOT NULL,
    pocket_volume_a3    REAL    NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_structures_gene ON protein_structures(tax_id, gene_symbol);
";

const SQL_CREATE_DATASET_MANIFESTS: &str = "
CREATE TABLE IF NOT EXISTS dataset_manifests (
    manifest_id   TEXT    PRIMARY KEY,
    tax_id        INTEGER NOT NULL REFERENCES organisms(tax_id),
    assembly_name TEXT    NOT NULL,
    provider      TEXT    NOT NULL,
    source_sha256 TEXT    NOT NULL,
    ingested_at   TEXT    NOT NULL
);
";

const SQL_CREATE_REGULATORY_EDGES: &str = "
CREATE TABLE IF NOT EXISTS regulatory_edges (
    edge_id         TEXT    PRIMARY KEY,
    tax_id          INTEGER NOT NULL REFERENCES organisms(tax_id),
    tf_symbol       TEXT    NOT NULL,
    target_symbol   TEXT    NOT NULL,
    interaction     TEXT    NOT NULL,
    weight          REAL    NOT NULL,
    evidence_source TEXT    NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_reg_edges_tax ON regulatory_edges(tax_id);
";

const SQL_CREATE_MULTIOMICS_TRACKS: &str = "
CREATE TABLE IF NOT EXISTS multiomics_tracks (
    tax_id                    INTEGER NOT NULL REFERENCES organisms(tax_id),
    assembly_name             TEXT    NOT NULL,
    sample_id                 TEXT    NOT NULL,
    gene_symbol               TEXT    NOT NULL,
    rna_tpm                   REAL    NOT NULL,
    atac_accessibility        REAL    NOT NULL,
    promoter_methylation_beta REAL    NOT NULL,
    h3k27ac_chip              REAL    NOT NULL,
    h3k27ac_input             REAL    NOT NULL,
    enhancer_distance_bp      INTEGER NOT NULL,
    PRIMARY KEY (tax_id, assembly_name, sample_id, gene_symbol)
);
";

const SQL_CREATE_SPATIAL_OMICS: &str = "
CREATE TABLE IF NOT EXISTS spatial_cells (
    dataset_id TEXT NOT NULL,
    cell_id    TEXT NOT NULL,
    x_um       REAL NOT NULL,
    y_um       REAL NOT NULL,
    PRIMARY KEY (dataset_id, cell_id)
);
CREATE TABLE IF NOT EXISTS spatial_transcripts (
    dataset_id  TEXT NOT NULL,
    spot_id     TEXT NOT NULL,
    gene_symbol TEXT NOT NULL,
    x_um        REAL NOT NULL,
    y_um        REAL NOT NULL,
    PRIMARY KEY (dataset_id, spot_id)
);
CREATE INDEX IF NOT EXISTS idx_spatial_transcripts_ds ON spatial_transcripts(dataset_id);
";

const SQL_CREATE_SINGLECELL_PHARMA: &str = "
CREATE TABLE IF NOT EXISTS single_cell_atlas (
    dataset_id     TEXT NOT NULL,
    cell_id        TEXT NOT NULL,
    cell_type      TEXT NOT NULL,
    expression_csv TEXT NOT NULL,
    PRIMARY KEY (dataset_id, cell_id)
);
CREATE TABLE IF NOT EXISTS pharma_compounds (
    compound_id      TEXT PRIMARY KEY,
    compound_name    TEXT NOT NULL,
    target_gene      TEXT NOT NULL,
    ec50_um          REAL NOT NULL,
    td50_um          REAL NOT NULL,
    hill_coefficient REAL NOT NULL,
    half_life_hours  REAL NOT NULL,
    bioavailability  REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS synthetic_gametes (
    gamete_id           TEXT PRIMARY KEY,
    tax_id              INTEGER NOT NULL REFERENCES organisms(tax_id),
    gamete_type         TEXT NOT NULL,
    follicle_stage      TEXT NOT NULL,
    diameter_um         REAL NOT NULL,
    mitochondrial_count INTEGER NOT NULL,
    imprint_fidelity    REAL NOT NULL,
    aneuploidy_risk     REAL NOT NULL,
    overall_score       REAL NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_synthetic_gametes_tax ON synthetic_gametes(tax_id);
";

const SQL_CREATE_LONGEVITY_TARGETS: &str = "
CREATE TABLE IF NOT EXISTS longevity_targets (
    symbol           TEXT    NOT NULL,
    tax_id           INTEGER NOT NULL,
    priority         INTEGER NOT NULL,
    category         TEXT    NOT NULL,
    default_assembly TEXT,
    notes            TEXT,
    PRIMARY KEY (tax_id, symbol)
);
CREATE INDEX IF NOT EXISTS idx_longevity_tax ON longevity_targets(tax_id, priority);
CREATE INDEX IF NOT EXISTS idx_longevity_symbol ON longevity_targets(symbol);
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
            "{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}",
            SQL_CREATE_ORGANISMS,
            SQL_CREATE_REFERENCE_GENOMES,
            SQL_CREATE_CHROMOSOMES,
            SQL_CREATE_ASSEMBLIES,
            SQL_CREATE_SEQUENCES,
            SQL_CREATE_GENES,
            SQL_CREATE_TRANSCRIPTS,
            SQL_CREATE_EXONS,
            SQL_CREATE_PROTEIN_ENTRIES,
            SQL_CREATE_PROTEIN_DOMAINS,
            SQL_CREATE_VARIANT_RECORDS,
            SQL_CREATE_CLINICAL_VARIANTS,
            SQL_CREATE_EPIGENETIC_FEATURES,
            SQL_CREATE_PROTEIN_STRUCTURES,
            SQL_CREATE_DATASET_MANIFESTS,
            SQL_CREATE_REGULATORY_EDGES,
            SQL_CREATE_MULTIOMICS_TRACKS,
            SQL_CREATE_SPATIAL_OMICS,
            SQL_CREATE_SINGLECELL_PHARMA,
            SQL_CREATE_LONGEVITY_TARGETS,
        );
        conn.execute_batch(&schema_batch)?;
        Self::seed_default_longevity_targets_conn(&conn)?;
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
    // Assemblies & Sequences (Section 4.1)
    // -----------------------------------------------------------------------

    pub fn insert_assembly(&self, assembly: &Assembly) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO assemblies (tax_id, assembly_name, source, created_at, notes)
             VALUES (?1, ?2, ?3, ?4, ?5)
             ON CONFLICT (tax_id, assembly_name) DO UPDATE SET
                source = COALESCE(excluded.source, assemblies.source),
                created_at = COALESCE(excluded.created_at, assemblies.created_at),
                notes = COALESCE(excluded.notes, assemblies.notes)",
            params![
                assembly.tax_id,
                assembly.assembly_name,
                assembly.source,
                assembly.created_at,
                assembly.notes,
            ],
        )?;
        Ok(())
    }

    pub fn get_assembly(&self, tax_id: u32, assembly_name: &str) -> Result<Option<Assembly>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, assembly_name, source, created_at, notes
             FROM assemblies
             WHERE tax_id = ?1 AND assembly_name = ?2",
        )?;
        let asm = stmt
            .query_row(params![tax_id, assembly_name], |row| {
                Ok(Assembly {
                    tax_id: row.get(0)?,
                    assembly_name: row.get(1)?,
                    source: row.get(2)?,
                    created_at: row.get(3)?,
                    notes: row.get(4)?,
                })
            })
            .optional()?;
        Ok(asm)
    }

    pub fn list_assemblies(&self, tax_id: u32) -> Result<Vec<Assembly>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, assembly_name, source, created_at, notes
             FROM assemblies
             WHERE tax_id = ?1
             ORDER BY assembly_name",
        )?;
        let rows = stmt.query_map(params![tax_id], |row| {
            Ok(Assembly {
                tax_id: row.get(0)?,
                assembly_name: row.get(1)?,
                source: row.get(2)?,
                created_at: row.get(3)?,
                notes: row.get(4)?,
            })
        })?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    pub fn insert_sequence(&self, seq: &SequenceRef) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO sequences (tax_id, assembly_name, seq_id, length_bp, path, sha256, source_uri)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)
             ON CONFLICT (tax_id, assembly_name, seq_id) DO UPDATE SET
                length_bp = excluded.length_bp,
                path = excluded.path,
                sha256 = excluded.sha256,
                source_uri = excluded.source_uri",
            params![
                seq.tax_id,
                seq.assembly_name,
                seq.seq_id,
                seq.length_bp,
                seq.path,
                seq.sha256,
                seq.source_uri,
            ],
        )?;
        Ok(())
    }

    pub fn get_sequence_ref(
        &self,
        tax_id: u32,
        assembly_name: &str,
        seq_id: &str,
    ) -> Result<Option<SequenceRef>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, assembly_name, seq_id, length_bp, path, sha256, source_uri
             FROM sequences
             WHERE tax_id = ?1 AND assembly_name = ?2 AND seq_id = ?3",
        )?;
        let seq = stmt
            .query_row(params![tax_id, assembly_name, seq_id], |row| {
                Ok(SequenceRef {
                    tax_id: row.get(0)?,
                    assembly_name: row.get(1)?,
                    seq_id: row.get(2)?,
                    length_bp: row.get(3)?,
                    path: row.get(4)?,
                    sha256: row.get(5)?,
                    source_uri: row.get(6)?,
                })
            })
            .optional()?;
        Ok(seq)
    }

    pub fn list_sequences(&self, tax_id: u32, assembly_name: &str) -> Result<Vec<SequenceRef>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, assembly_name, seq_id, length_bp, path, sha256, source_uri
             FROM sequences
             WHERE tax_id = ?1 AND assembly_name = ?2
             ORDER BY seq_id",
        )?;
        let rows = stmt.query_map(params![tax_id, assembly_name], |row| {
            Ok(SequenceRef {
                tax_id: row.get(0)?,
                assembly_name: row.get(1)?,
                seq_id: row.get(2)?,
                length_bp: row.get(3)?,
                path: row.get(4)?,
                sha256: row.get(5)?,
                source_uri: row.get(6)?,
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

    /// Inserts a slice of [`Gene`] records within a single atomic SQLite transaction.
    pub fn insert_genes_batch(&self, genes: &[Gene]) -> Result<()> {
        if genes.is_empty() {
            return Ok(());
        }
        let mut conn = self.conn.lock().unwrap();
        let tx = conn.transaction()?;
        {
            let mut stmt = tx.prepare(
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
            )?;
            for gene in genes {
                let biotype_str = match &gene.biotype {
                    Biotype::ProteinCoding => "protein_coding",
                    Biotype::LongNonCodingRna => "lncRNA",
                    Biotype::MicroRna => "miRNA",
                    Biotype::SmallNuclearRna => "snRNA",
                    Biotype::Pseudogene => "pseudogene",
                    Biotype::RibosomalRna => "rRNA",
                    Biotype::Other(s) => s.as_str(),
                };
                let strand_str = match gene.strand {
                    Strand::Forward => "+",
                    Strand::Reverse => "-",
                };
                stmt.execute(params![
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
                ])?;
            }
        }
        tx.commit()?;
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

    /// Returns all genes on chromosome `chr` scoped to `tax_id` and `assembly`, ordered by start coordinate.
    pub fn find_genes_by_chromosome(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
    ) -> Result<Vec<Gene>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id
             FROM genes
             WHERE tax_id = ?1 AND assembly_name = ?2 AND chr = ?3
             ORDER BY start",
        )?;
        let rows = stmt.query_map(params![tax_id, assembly, chr], Self::row_to_gene)?;
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

    fn row_to_gene_feature(row: &rusqlite::Row) -> rusqlite::Result<GeneFeature> {
        let strand_raw: String = row.get(7)?;
        let strand = if strand_raw == "-" {
            Strand::Reverse
        } else {
            Strand::Forward
        };
        Ok(GeneFeature {
            tax_id: row.get(0)?,
            assembly_name: row.get(1)?,
            gene_id: row.get(2)?,
            symbol: row.get(3)?,
            chrom: row.get(4)?,
            start_bp: row.get(5)?,
            end_bp: row.get(6)?,
            strand,
            feature_type: row.get(8)?,
            source: row.get(9)?,
        })
    }

    /// Queries a gene feature by symbol (or gene_id) scoped by `tax_id` and `assembly` (Section 5 DB4).
    pub fn get_gene_feature_by_symbol(
        &self,
        tax_id: u32,
        assembly: &str,
        symbol: &str,
    ) -> Result<Option<GeneFeature>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, assembly_name, gene_id, symbol, chr, start, end, strand, feature_type, source
             FROM genes
             WHERE tax_id = ?1 AND assembly_name = ?2 AND (symbol = ?3 OR gene_id = ?3)
             ORDER BY CASE WHEN feature_type = 'gene' THEN 0 ELSE 1 END
             LIMIT 1",
        )?;
        let feat = stmt
            .query_row(params![tax_id, assembly, symbol], Self::row_to_gene_feature)
            .optional()?;
        Ok(feat)
    }

    /// Lists all gene features overlapping interval `[start, end)` on chromosome `chrom` (Section 5 DB4).
    pub fn list_gene_features_in_region(
        &self,
        tax_id: u32,
        assembly: &str,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<GeneFeature>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, assembly_name, gene_id, symbol, chr, start, end, strand, feature_type, source
             FROM genes
             WHERE tax_id = ?1 AND assembly_name = ?2 AND chr = ?3 AND start < ?5 AND end > ?4
             ORDER BY start",
        )?;
        let rows = stmt.query_map(
            params![tax_id, assembly, chrom, start, end],
            Self::row_to_gene_feature,
        )?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    /// Inserts a [`GeneFeature`] into the `genes` table.
    pub fn insert_gene_feature(&self, feat: &GeneFeature) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let strand_str = match feat.strand {
            Strand::Forward => "+",
            Strand::Reverse => "-",
        };
        conn.execute(
            "INSERT INTO genes (gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id, feature_type, source)
             VALUES (?1, ?2, ?3, ?4, 'protein_coding', ?5, ?6, ?7, ?8, ?9, ?10, ?11)
             ON CONFLICT (gene_id, assembly_name) DO UPDATE SET
                symbol = excluded.symbol,
                name = excluded.name,
                chr = excluded.chr,
                start = excluded.start,
                end = excluded.end,
                strand = excluded.strand,
                tax_id = excluded.tax_id,
                feature_type = excluded.feature_type,
                source = excluded.source",
            params![
                feat.gene_id,
                feat.assembly_name,
                feat.symbol,
                feat.symbol,
                feat.chrom,
                feat.start_bp,
                feat.end_bp,
                strand_str,
                feat.tax_id,
                feat.feature_type,
                feat.source,
            ],
        )?;
        Ok(())
    }

    /// Batch-inserts [`GeneFeature`] records within an atomic transaction.
    pub fn insert_gene_features_batch(&self, feats: &[GeneFeature]) -> Result<()> {
        let mut conn = self.conn.lock().unwrap();
        let tx = conn.transaction()?;
        {
            let mut stmt = tx.prepare(
                "INSERT INTO genes (gene_id, assembly_name, symbol, name, biotype, chr, start, end, strand, tax_id, feature_type, source)
                 VALUES (?1, ?2, ?3, ?4, 'protein_coding', ?5, ?6, ?7, ?8, ?9, ?10, ?11)
                 ON CONFLICT (gene_id, assembly_name) DO UPDATE SET
                    symbol = excluded.symbol,
                    name = excluded.name,
                    chr = excluded.chr,
                    start = excluded.start,
                    end = excluded.end,
                    strand = excluded.strand,
                    tax_id = excluded.tax_id,
                    feature_type = excluded.feature_type,
                    source = excluded.source",
            )?;
            for feat in feats {
                let strand_str = match feat.strand {
                    Strand::Forward => "+",
                    Strand::Reverse => "-",
                };
                stmt.execute(params![
                    feat.gene_id,
                    feat.assembly_name,
                    feat.symbol,
                    feat.symbol,
                    feat.chrom,
                    feat.start_bp,
                    feat.end_bp,
                    strand_str,
                    feat.tax_id,
                    feat.feature_type,
                    feat.source,
                ])?;
            }
        }
        tx.commit()?;
        Ok(())
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

    /// Inserts a slice of [`Transcript`] records within a single atomic SQLite transaction.
    pub fn insert_transcripts_batch(&self, transcripts: &[Transcript]) -> Result<()> {
        if transcripts.is_empty() {
            return Ok(());
        }
        let mut conn = self.conn.lock().unwrap();
        let tx = conn.transaction()?;
        {
            let mut stmt = tx.prepare(
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
            )?;
            for transcript in transcripts {
                let biotype_str = match &transcript.biotype {
                    TranscriptBiotype::MRNA => "mRNA",
                    TranscriptBiotype::NCRNA => "ncRNA",
                    TranscriptBiotype::Pseudogene => "pseudogene",
                    TranscriptBiotype::Other(s) => s.as_str(),
                };
                let strand_str = match transcript.strand {
                    Strand::Forward => "+",
                    Strand::Reverse => "-",
                };
                stmt.execute(params![
                    transcript.transcript_id,
                    transcript.gene_id,
                    transcript.assembly_name,
                    biotype_str,
                    transcript.chr,
                    transcript.start,
                    transcript.end,
                    strand_str,
                    transcript.exon_count
                ])?;
            }
        }
        tx.commit()?;
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

    /// Inserts a slice of [`Exon`] records within a single atomic SQLite transaction.
    pub fn insert_exons_batch(&self, exons: &[Exon]) -> Result<()> {
        if exons.is_empty() {
            return Ok(());
        }
        let mut conn = self.conn.lock().unwrap();
        let tx = conn.transaction()?;
        {
            let mut stmt = tx.prepare(
                "INSERT INTO exons (exon_id, transcript_id, assembly_name, chr, start, end, strand, exon_number)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)
                 ON CONFLICT (exon_id, transcript_id, assembly_name) DO UPDATE SET
                    chr = excluded.chr,
                    start = excluded.start,
                    end = excluded.end,
                    strand = excluded.strand,
                    exon_number = excluded.exon_number",
            )?;
            for exon in exons {
                let strand_str = match exon.strand {
                    Strand::Forward => "+",
                    Strand::Reverse => "-",
                };
                stmt.execute(params![
                    exon.exon_id,
                    exon.transcript_id,
                    exon.assembly_name,
                    exon.chr,
                    exon.start,
                    exon.end,
                    strand_str,
                    exon.exon_number
                ])?;
            }
        }
        tx.commit()?;
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

    /// Returns all transcripts on chromosome `chr` within `assembly`, ordered by start coordinate.
    pub fn find_transcripts_by_chromosome(
        &self,
        assembly: &str,
        chr: &str,
    ) -> Result<Vec<Transcript>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT transcript_id, gene_id, assembly_name, biotype, chr, start, end, strand, exon_count
             FROM transcripts WHERE assembly_name = ?1 AND chr = ?2 ORDER BY start",
        )?;
        let rows = stmt.query_map(params![assembly, chr], Self::row_to_transcript)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    /// Tunes SQLite PRAGMAs for high-throughput bulk insertion (larger cache and memory temp store).
    pub fn tune_for_bulk_loading(&self) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute_batch(
            "PRAGMA synchronous = NORMAL;
             PRAGMA cache_size = -64000;
             PRAGMA temp_store = MEMORY;",
        )?;
        Ok(())
    }

    /// Restores standard SQLite operational PRAGMAs.
    pub fn restore_normal_pragmas(&self) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute_batch(
            "PRAGMA synchronous = NORMAL;
             PRAGMA cache_size = -2000;
             PRAGMA temp_store = DEFAULT;",
        )?;
        Ok(())
    }

    /// Executes closure `f` while SQLite PRAGMAs are tuned for high-throughput bulk loading.
    pub fn with_bulk_mode<F, T>(&self, f: F) -> Result<T>
    where
        F: FnOnce() -> Result<T>,
    {
        self.tune_for_bulk_loading()?;
        let res = f();
        let _ = self.restore_normal_pragmas();
        res
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

    /// Inserts a slice of [`VariantRecord`] records within a single atomic SQLite transaction.
    pub fn insert_variants_batch(&self, variants: &[VariantRecord]) -> Result<()> {
        if variants.is_empty() {
            return Ok(());
        }
        let mut conn = self.conn.lock().unwrap();
        let tx = conn.transaction()?;
        {
            let mut stmt = tx.prepare(
                "INSERT INTO variant_records (variant_id, assembly_name, chr, position, reference, alternate, variant_type, tax_id)
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)
                 ON CONFLICT (variant_id, assembly_name) DO UPDATE SET
                    chr = excluded.chr,
                    position = excluded.position,
                    reference = excluded.reference,
                    alternate = excluded.alternate,
                    variant_type = excluded.variant_type,
                    tax_id = excluded.tax_id",
            )?;
            for variant in variants {
                let type_str = match variant.variant_type {
                    VariantType::Snp => "SNP",
                    VariantType::Insertion => "Insertion",
                    VariantType::Deletion => "Deletion",
                    VariantType::Mnp => "MNP",
                    VariantType::StructuralVariant => "SV",
                };
                stmt.execute(params![
                    variant.variant_id,
                    variant.assembly_name,
                    variant.chr,
                    variant.position,
                    variant.reference,
                    variant.alternate,
                    type_str,
                    variant.tax_id
                ])?;
            }
        }
        tx.commit()?;
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

    // -----------------------------------------------------------------------
    // Epigenomics, 3D Protein Structures & Provenance Manifests (v0.4/v0.5)
    // -----------------------------------------------------------------------

    /// Inserts or updates an [`EpigeneticFeature`] record.
    pub fn insert_epigenetic_feature(&self, feat: &EpigeneticFeature) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let track_str = match feat.track_type {
            EpigeneticTrackType::CpgIsland => "CpgIsland",
            EpigeneticTrackType::HistoneMark => "HistoneMark",
            EpigeneticTrackType::ChromatinPeak => "ChromatinPeak",
            EpigeneticTrackType::RegulatoryElement => "RegulatoryElement",
        };
        conn.execute(
            "INSERT INTO epigenetic_features (feature_id, tax_id, assembly_name, chr, start, end, track_type, mark_label, signal_value)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)
             ON CONFLICT (feature_id, assembly_name) DO UPDATE SET
                chr = excluded.chr,
                start = excluded.start,
                end = excluded.end,
                track_type = excluded.track_type,
                mark_label = excluded.mark_label,
                signal_value = excluded.signal_value",
            params![
                feat.feature_id,
                feat.tax_id,
                feat.assembly_name,
                feat.chr,
                feat.start,
                feat.end,
                track_str,
                feat.mark_label,
                feat.signal_value
            ],
        )?;
        Ok(())
    }

    /// Queries all [`EpigeneticFeature`] records overlapping `[start, end)` on `chr`.
    pub fn find_epigenetic_features_by_region(
        &self,
        tax_id: u32,
        assembly_name: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<EpigeneticFeature>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT feature_id, tax_id, assembly_name, chr, start, end, track_type, mark_label, signal_value
             FROM epigenetic_features
             WHERE tax_id = ?1 AND assembly_name = ?2 AND chr = ?3 AND start < ?5 AND end > ?4
             ORDER BY start",
        )?;
        let rows = stmt.query_map(params![tax_id, assembly_name, chr, start, end], |row| {
            let track_str: String = row.get(6)?;
            let track_type = match track_str.as_str() {
                "CpgIsland" => EpigeneticTrackType::CpgIsland,
                "HistoneMark" => EpigeneticTrackType::HistoneMark,
                "ChromatinPeak" => EpigeneticTrackType::ChromatinPeak,
                _ => EpigeneticTrackType::RegulatoryElement,
            };
            Ok(EpigeneticFeature {
                feature_id: row.get(0)?,
                tax_id: row.get(1)?,
                assembly_name: row.get(2)?,
                chr: row.get(3)?,
                start: row.get(4)?,
                end: row.get(5)?,
                track_type,
                mark_label: row.get(7)?,
                signal_value: row.get(8)?,
            })
        })?;
        let mut out = Vec::new();
        for r in rows {
            out.push(r?);
        }
        Ok(out)
    }

    /// Inserts or updates a 3D [`ProteinStructureRecord`].
    pub fn insert_protein_structure(&self, rec: &ProteinStructureRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO protein_structures (structure_id, gene_symbol, tax_id, resolution_angstrom, residue_count, pocket_volume_a3)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)
             ON CONFLICT (structure_id) DO UPDATE SET
                gene_symbol = excluded.gene_symbol,
                resolution_angstrom = excluded.resolution_angstrom,
                residue_count = excluded.residue_count,
                pocket_volume_a3 = excluded.pocket_volume_a3",
            params![
                rec.structure_id,
                rec.gene_symbol,
                rec.tax_id,
                rec.resolution_angstrom,
                rec.residue_count,
                rec.pocket_volume_a3
            ],
        )?;
        Ok(())
    }

    /// Finds 3D [`ProteinStructureRecord`] entries for `(tax_id, gene_symbol)`.
    pub fn find_protein_structures_by_gene(
        &self,
        tax_id: u32,
        gene_symbol: &str,
    ) -> Result<Vec<ProteinStructureRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT structure_id, gene_symbol, tax_id, resolution_angstrom, residue_count, pocket_volume_a3
             FROM protein_structures WHERE tax_id = ?1 AND gene_symbol = ?2 ORDER BY resolution_angstrom",
        )?;
        let rows = stmt.query_map(params![tax_id, gene_symbol], |row| {
            Ok(ProteinStructureRecord {
                structure_id: row.get(0)?,
                gene_symbol: row.get(1)?,
                tax_id: row.get(2)?,
                resolution_angstrom: row.get(3)?,
                residue_count: row.get(4)?,
                pocket_volume_a3: row.get(5)?,
            })
        })?;
        let mut out = Vec::new();
        for r in rows {
            out.push(r?);
        }
        Ok(out)
    }

    /// Records a cryptographic [`DatasetManifest`] provenance entry.
    pub fn insert_dataset_manifest(&self, manifest: &DatasetManifest) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO dataset_manifests (manifest_id, tax_id, assembly_name, provider, source_sha256, ingested_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)
             ON CONFLICT (manifest_id) DO UPDATE SET
                source_sha256 = excluded.source_sha256,
                ingested_at = excluded.ingested_at",
            params![
                manifest.manifest_id,
                manifest.tax_id,
                manifest.assembly_name,
                manifest.provider,
                manifest.source_sha256,
                manifest.ingested_at
            ],
        )?;
        Ok(())
    }

    /// Retrieves a [`DatasetManifest`] by `manifest_id`.
    pub fn get_dataset_manifest(&self, manifest_id: &str) -> Result<Option<DatasetManifest>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT manifest_id, tax_id, assembly_name, provider, source_sha256, ingested_at
             FROM dataset_manifests WHERE manifest_id = ?1",
        )?;
        stmt.query_row(params![manifest_id], |row| {
            Ok(DatasetManifest {
                manifest_id: row.get(0)?,
                tax_id: row.get(1)?,
                assembly_name: row.get(2)?,
                provider: row.get(3)?,
                source_sha256: row.get(4)?,
                ingested_at: row.get(5)?,
            })
        })
        .optional()
        .map_err(Into::into)
    }

    // -----------------------------------------------------------------------
    // Regulatory Networks, Multi-Omics Tracks & Spatial Omics (v0.4 / th v0.8)
    // -----------------------------------------------------------------------

    /// Inserts or updates a [`RegulatoryEdgeRecord`] in `regulatory_edges`.
    pub fn insert_regulatory_edge(&self, edge: &RegulatoryEdgeRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let mode_str = match edge.interaction {
            scies_bio_th::regulatory::RegulationType::Activation => "Activation",
            scies_bio_th::regulatory::RegulationType::Repression => "Repression",
        };
        conn.execute(
            "INSERT INTO regulatory_edges (edge_id, tax_id, tf_symbol, target_symbol, interaction, weight, evidence_source)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)
             ON CONFLICT (edge_id) DO UPDATE SET
                interaction = excluded.interaction,
                weight = excluded.weight,
                evidence_source = excluded.evidence_source",
            params![
                edge.edge_id,
                edge.tax_id,
                edge.tf_symbol,
                edge.target_symbol,
                mode_str,
                edge.weight,
                edge.evidence_source
            ],
        )?;
        Ok(())
    }

    /// Lists all [`RegulatoryEdgeRecord`] interactions for `tax_id`.
    pub fn list_regulatory_edges(&self, tax_id: u32) -> Result<Vec<RegulatoryEdgeRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT edge_id, tax_id, tf_symbol, target_symbol, interaction, weight, evidence_source
             FROM regulatory_edges WHERE tax_id = ?1 ORDER BY edge_id",
        )?;
        let rows = stmt.query_map(params![tax_id], |row| {
            let mode_str: String = row.get(4)?;
            let interaction = if mode_str == "Repression" {
                scies_bio_th::regulatory::RegulationType::Repression
            } else {
                scies_bio_th::regulatory::RegulationType::Activation
            };
            Ok(RegulatoryEdgeRecord {
                edge_id: row.get(0)?,
                tax_id: row.get(1)?,
                tf_symbol: row.get(2)?,
                target_symbol: row.get(3)?,
                interaction,
                weight: row.get(5)?,
                evidence_source: row.get(6)?,
            })
        })?;
        let mut out = Vec::new();
        for r in rows {
            out.push(r?);
        }
        Ok(out)
    }

    /// Inserts or updates a [`MultiOmicsTrackRecord`] in `multiomics_tracks`.
    pub fn insert_multiomics_track(&self, rec: &MultiOmicsTrackRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO multiomics_tracks (tax_id, assembly_name, sample_id, gene_symbol, rna_tpm, atac_accessibility, promoter_methylation_beta, h3k27ac_chip, h3k27ac_input, enhancer_distance_bp)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)
             ON CONFLICT (tax_id, assembly_name, sample_id, gene_symbol) DO UPDATE SET
                rna_tpm = excluded.rna_tpm,
                atac_accessibility = excluded.atac_accessibility,
                promoter_methylation_beta = excluded.promoter_methylation_beta,
                h3k27ac_chip = excluded.h3k27ac_chip,
                h3k27ac_input = excluded.h3k27ac_input,
                enhancer_distance_bp = excluded.enhancer_distance_bp",
            params![
                rec.tax_id,
                rec.assembly_name,
                rec.sample_id,
                rec.gene_symbol,
                rec.rna_tpm,
                rec.atac_accessibility,
                rec.promoter_methylation_beta,
                rec.h3k27ac_chip,
                rec.h3k27ac_input,
                rec.enhancer_distance_bp
            ],
        )?;
        Ok(())
    }

    /// Fetches a [`MultiOmicsTrackRecord`] for `(tax_id, assembly_name, sample_id, gene_symbol)`.
    pub fn get_multiomics_track(
        &self,
        tax_id: u32,
        assembly_name: &str,
        sample_id: &str,
        gene_symbol: &str,
    ) -> Result<Option<MultiOmicsTrackRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT tax_id, assembly_name, sample_id, gene_symbol, rna_tpm, atac_accessibility, promoter_methylation_beta, h3k27ac_chip, h3k27ac_input, enhancer_distance_bp
             FROM multiomics_tracks
             WHERE tax_id = ?1 AND assembly_name = ?2 AND sample_id = ?3 AND gene_symbol = ?4",
        )?;
        stmt.query_row(
            params![tax_id, assembly_name, sample_id, gene_symbol],
            |row| {
                Ok(MultiOmicsTrackRecord {
                    tax_id: row.get(0)?,
                    assembly_name: row.get(1)?,
                    sample_id: row.get(2)?,
                    gene_symbol: row.get(3)?,
                    rna_tpm: row.get(4)?,
                    atac_accessibility: row.get(5)?,
                    promoter_methylation_beta: row.get(6)?,
                    h3k27ac_chip: row.get(7)?,
                    h3k27ac_input: row.get(8)?,
                    enhancer_distance_bp: row.get(9)?,
                })
            },
        )
        .optional()
        .map_err(Into::into)
    }

    /// Inserts or updates a [`SpatialCellRecord`] in `spatial_cells`.
    pub fn insert_spatial_cell(&self, cell: &SpatialCellRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO spatial_cells (dataset_id, cell_id, x_um, y_um)
             VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT (dataset_id, cell_id) DO UPDATE SET
                x_um = excluded.x_um,
                y_um = excluded.y_um",
            params![cell.dataset_id, cell.cell_id, cell.x_um, cell.y_um],
        )?;
        Ok(())
    }

    /// Lists all [`SpatialCellRecord`] entries for `dataset_id`.
    pub fn list_spatial_cells(&self, dataset_id: &str) -> Result<Vec<SpatialCellRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT dataset_id, cell_id, x_um, y_um FROM spatial_cells WHERE dataset_id = ?1 ORDER BY cell_id",
        )?;
        let rows = stmt.query_map(params![dataset_id], |row| {
            Ok(SpatialCellRecord {
                dataset_id: row.get(0)?,
                cell_id: row.get(1)?,
                x_um: row.get(2)?,
                y_um: row.get(3)?,
            })
        })?;
        let mut out = Vec::new();
        for r in rows {
            out.push(r?);
        }
        Ok(out)
    }

    /// Inserts or updates a [`SpatialTranscriptRecord`] in `spatial_transcripts`.
    pub fn insert_spatial_transcript(&self, spot: &SpatialTranscriptRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO spatial_transcripts (dataset_id, spot_id, gene_symbol, x_um, y_um)
             VALUES (?1, ?2, ?3, ?4, ?5)
             ON CONFLICT (dataset_id, spot_id) DO UPDATE SET
                gene_symbol = excluded.gene_symbol,
                x_um = excluded.x_um,
                y_um = excluded.y_um",
            params![
                spot.dataset_id,
                spot.spot_id,
                spot.gene_symbol,
                spot.x_um,
                spot.y_um
            ],
        )?;
        Ok(())
    }

    /// Lists all [`SpatialTranscriptRecord`] spots for `dataset_id`.
    pub fn list_spatial_transcripts(
        &self,
        dataset_id: &str,
    ) -> Result<Vec<SpatialTranscriptRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT dataset_id, spot_id, gene_symbol, x_um, y_um FROM spatial_transcripts WHERE dataset_id = ?1 ORDER BY spot_id",
        )?;
        let rows = stmt.query_map(params![dataset_id], |row| {
            Ok(SpatialTranscriptRecord {
                dataset_id: row.get(0)?,
                spot_id: row.get(1)?,
                gene_symbol: row.get(2)?,
                x_um: row.get(3)?,
                y_um: row.get(4)?,
            })
        })?;
        let mut out = Vec::new();
        for r in rows {
            out.push(r?);
        }
        Ok(out)
    }

    // -----------------------------------------------------------------------
    // Single-Cell Atlas & Pharmacokinetic Compounds (v0.5 / th v0.9)
    // -----------------------------------------------------------------------

    /// Inserts or updates a [`SingleCellAtlasRecord`] in `single_cell_atlas`.
    pub fn insert_single_cell_record(&self, rec: &SingleCellAtlasRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        let csv = rec
            .expression
            .iter()
            .map(|v| format!("{:.6}", v))
            .collect::<Vec<_>>()
            .join(",");
        conn.execute(
            "INSERT INTO single_cell_atlas (dataset_id, cell_id, cell_type, expression_csv)
             VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT (dataset_id, cell_id) DO UPDATE SET
                cell_type = excluded.cell_type,
                expression_csv = excluded.expression_csv",
            params![rec.dataset_id, rec.cell_id, rec.cell_type, csv],
        )?;
        Ok(())
    }

    /// Lists all [`SingleCellAtlasRecord`] entries for `dataset_id`.
    pub fn list_single_cell_records(
        &self,
        dataset_id: &str,
    ) -> Result<Vec<SingleCellAtlasRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT dataset_id, cell_id, cell_type, expression_csv
             FROM single_cell_atlas WHERE dataset_id = ?1 ORDER BY cell_id",
        )?;
        let rows = stmt.query_map(params![dataset_id], |row| {
            let csv: String = row.get(3)?;
            let expression: Vec<f64> = csv
                .split(',')
                .filter(|s| !s.is_empty())
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            Ok(SingleCellAtlasRecord {
                dataset_id: row.get(0)?,
                cell_id: row.get(1)?,
                cell_type: row.get(2)?,
                expression,
            })
        })?;
        let mut out = Vec::new();
        for r in rows {
            out.push(r?);
        }
        Ok(out)
    }

    /// Inserts or updates a [`PharmaCompoundRecord`] in `pharma_compounds`.
    pub fn insert_pharma_compound(&self, rec: &PharmaCompoundRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO pharma_compounds (compound_id, compound_name, target_gene, ec50_um, td50_um, hill_coefficient, half_life_hours, bioavailability)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)
             ON CONFLICT (compound_id) DO UPDATE SET
                compound_name = excluded.compound_name,
                target_gene = excluded.target_gene,
                ec50_um = excluded.ec50_um,
                td50_um = excluded.td50_um,
                hill_coefficient = excluded.hill_coefficient,
                half_life_hours = excluded.half_life_hours,
                bioavailability = excluded.bioavailability",
            params![
                rec.compound_id,
                rec.compound_name,
                rec.target_gene,
                rec.ec50_um,
                rec.td50_um,
                rec.hill_coefficient,
                rec.half_life_hours,
                rec.bioavailability
            ],
        )?;
        Ok(())
    }

    /// Retrieves a [`PharmaCompoundRecord`] by `compound_id`.
    pub fn get_pharma_compound(&self, compound_id: &str) -> Result<Option<PharmaCompoundRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT compound_id, compound_name, target_gene, ec50_um, td50_um, hill_coefficient, half_life_hours, bioavailability
             FROM pharma_compounds WHERE compound_id = ?1",
        )?;
        stmt.query_row(params![compound_id], |row| {
            Ok(PharmaCompoundRecord {
                compound_id: row.get(0)?,
                compound_name: row.get(1)?,
                target_gene: row.get(2)?,
                ec50_um: row.get(3)?,
                td50_um: row.get(4)?,
                hill_coefficient: row.get(5)?,
                half_life_hours: row.get(6)?,
                bioavailability: row.get(7)?,
            })
        })
        .optional()
        .map_err(Into::into)
    }

    /// Inserts or replaces a [`SyntheticGameteRecord`] in `synthetic_gametes`.
    pub fn insert_synthetic_gamete(&self, record: &SyntheticGameteRecord) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO synthetic_gametes (gamete_id, tax_id, gamete_type, follicle_stage, diameter_um, mitochondrial_count, imprint_fidelity, aneuploidy_risk, overall_score)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)
             ON CONFLICT (gamete_id) DO UPDATE SET
                tax_id = excluded.tax_id,
                gamete_type = excluded.gamete_type,
                follicle_stage = excluded.follicle_stage,
                diameter_um = excluded.diameter_um,
                mitochondrial_count = excluded.mitochondrial_count,
                imprint_fidelity = excluded.imprint_fidelity,
                aneuploidy_risk = excluded.aneuploidy_risk,
                overall_score = excluded.overall_score",
            params![
                record.gamete_id,
                record.tax_id,
                record.gamete_type,
                record.follicle_stage,
                record.diameter_um,
                record.mitochondrial_count as i64,
                record.imprint_fidelity,
                record.aneuploidy_risk,
                record.overall_score,
            ],
        )?;
        Ok(())
    }

    /// Retrieves a [`SyntheticGameteRecord`] by `gamete_id`.
    pub fn get_synthetic_gamete(&self, gamete_id: &str) -> Result<Option<SyntheticGameteRecord>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT gamete_id, tax_id, gamete_type, follicle_stage, diameter_um, mitochondrial_count, imprint_fidelity, aneuploidy_risk, overall_score
             FROM synthetic_gametes WHERE gamete_id = ?1",
        )?;
        stmt.query_row(params![gamete_id], |row| {
            let count_i64: i64 = row.get(5)?;
            Ok(SyntheticGameteRecord {
                gamete_id: row.get(0)?,
                tax_id: row.get(1)?,
                gamete_type: row.get(2)?,
                follicle_stage: row.get(3)?,
                diameter_um: row.get(4)?,
                mitochondrial_count: count_i64.max(0) as u64,
                imprint_fidelity: row.get(6)?,
                aneuploidy_risk: row.get(7)?,
                overall_score: row.get(8)?,
            })
        })
        .optional()
        .map_err(Into::into)
    }

    // -----------------------------------------------------------------------
    // Longevity Targets (Section 4.1 & Sprint DB-S3)
    // -----------------------------------------------------------------------

    fn row_to_longevity_target(row: &rusqlite::Row) -> rusqlite::Result<LongevityTarget> {
        Ok(LongevityTarget {
            symbol: row.get(0)?,
            tax_id: row.get(1)?,
            priority: row.get(2)?,
            category: row.get(3)?,
            default_assembly: row.get(4)?,
            notes: row.get(5)?,
        })
    }

    /// Inserts or updates a [`LongevityTarget`] in the catalog.
    pub fn insert_longevity_target(&self, target: &LongevityTarget) -> Result<()> {
        let conn = self.conn.lock().unwrap();
        conn.execute(
            "INSERT INTO longevity_targets (symbol, tax_id, priority, category, default_assembly, notes)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)
             ON CONFLICT (tax_id, symbol) DO UPDATE SET
                priority = excluded.priority,
                category = excluded.category,
                default_assembly = excluded.default_assembly,
                notes = excluded.notes",
            params![
                target.symbol,
                target.tax_id,
                target.priority,
                target.category,
                target.default_assembly,
                target.notes,
            ],
        )?;
        Ok(())
    }

    /// Retrieves a [`LongevityTarget`] for an organism `tax_id` and `symbol`.
    pub fn get_longevity_target(&self, tax_id: u32, symbol: &str) -> Result<Option<LongevityTarget>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT symbol, tax_id, priority, category, default_assembly, notes
             FROM longevity_targets
             WHERE tax_id = ?1 AND symbol = ?2",
        )?;
        let target = stmt
            .query_row(params![tax_id, symbol], Self::row_to_longevity_target)
            .optional()?;
        Ok(target)
    }

    /// Retrieves a [`LongevityTarget`] by `symbol` across all organisms (preferring human tax_id 9606).
    pub fn get_longevity_target_by_symbol(&self, symbol: &str) -> Result<Option<LongevityTarget>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT symbol, tax_id, priority, category, default_assembly, notes
             FROM longevity_targets
             WHERE symbol = ?1
             ORDER BY CASE WHEN tax_id = 9606 THEN 0 ELSE 1 END
             LIMIT 1",
        )?;
        let target = stmt
            .query_row(params![symbol], Self::row_to_longevity_target)
            .optional()?;
        Ok(target)
    }

    /// Lists all registered longevity targets for `tax_id` sorted by priority ASC, symbol ASC.
    pub fn list_longevity_targets(&self, tax_id: u32) -> Result<Vec<LongevityTarget>> {
        let conn = self.conn.lock().unwrap();
        let mut stmt = conn.prepare(
            "SELECT symbol, tax_id, priority, category, default_assembly, notes
             FROM longevity_targets
             WHERE tax_id = ?1
             ORDER BY priority ASC, symbol ASC",
        )?;
        let rows = stmt.query_map(params![tax_id], Self::row_to_longevity_target)?;
        let mut list = Vec::new();
        for r in rows {
            list.push(r?);
        }
        Ok(list)
    }

    /// Internal helper that seeds the default 10 essential human longevity targets (tax_id 9606) using a raw connection.
    pub(crate) fn seed_default_longevity_targets_conn(conn: &Connection) -> Result<usize> {
        let seed_targets: [(&str, u32, u32, &str, &str, &str); 10] = [
            ("TERT", 9606, 1, "telomere", "T2T-CHM13v2.0", "Telomerase reverse transcriptase; primary catalytic subunit for telomere maintenance and elongation"),
            ("TERC", 9606, 1, "telomere", "T2T-CHM13v2.0", "Telomerase RNA component; non-coding template RNA essential for telomerase activity"),
            ("SIRT6", 9606, 1, "dna_repair", "T2T-CHM13v2.0", "Sirtuin 6; NAD+-dependent deacetylase regulating DSB repair, chromatin stability, and glucose metabolism"),
            ("FOXO3", 9606, 2, "metabolic", "T2T-CHM13v2.0", "Forkhead box O3; transcription factor strongly linked to human exceptional longevity GWAS"),
            ("TP53", 9606, 1, "tumor_suppressor", "T2T-CHM13v2.0", "Tumor protein p53; genome guardian, safeguards cellular senescence and partial reprogramming safety"),
            ("LMNA", 9606, 2, "proteostasis", "T2T-CHM13v2.0", "Lamin A/C; nuclear lamina structural integrity, target in Hutchinson-Gilford progeria syndrome"),
            ("CDKN2A", 9606, 2, "tumor_suppressor", "T2T-CHM13v2.0", "Cyclin-dependent kinase inhibitor 2A (p16INK4a/p14ARF); biomarker and driver of senescent arrest"),
            ("ATM", 9606, 2, "dna_repair", "T2T-CHM13v2.0", "ATM serine/threonine kinase; master orchestrator of double-strand break response and DNA repair signaling"),
            ("PARP1", 9606, 2, "dna_repair", "T2T-CHM13v2.0", "Poly(ADP-ribose) polymerase 1; early DNA damage sensor and base excision repair regulator"),
            ("KLOTHO", 9606, 2, "metabolic", "T2T-CHM13v2.0", "Klotho alpha; anti-aging circulating endocrine factor suppressing IGF-1 and regulating phosphate homeostasis"),
        ];

        let mut inserted = 0;
        let mut stmt = conn.prepare(
            "INSERT OR IGNORE INTO longevity_targets (symbol, tax_id, priority, category, default_assembly, notes)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
        )?;
        for (symbol, tax_id, priority, category, assembly, notes) in seed_targets {
            let res = stmt.execute(params![symbol, tax_id, priority, category, assembly, notes])?;
            inserted += res;
        }
        Ok(inserted)
    }

    /// Seeds the default 10 essential human longevity targets (tax_id 9606).
    pub fn seed_default_longevity_targets(&self) -> Result<usize> {
        let conn = self.conn.lock().unwrap();
        Self::seed_default_longevity_targets_conn(&conn)
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

        // Synthetic Gamete CRUD check
        let gamete_rec = SyntheticGameteRecord {
            gamete_id: "OOCYTE_TEST_01".to_string(),
            tax_id: 9606,
            gamete_type: "Oocyte".to_string(),
            follicle_stage: "PreOvulatoryMII".to_string(),
            diameter_um: 118.5,
            mitochondrial_count: 220_000,
            imprint_fidelity: 0.985,
            aneuploidy_risk: 0.02,
            overall_score: 96.5,
        };
        store.insert_synthetic_gamete(&gamete_rec).unwrap();
        let fetched_gamete = store.get_synthetic_gamete("OOCYTE_TEST_01").unwrap().unwrap();
        assert_eq!(fetched_gamete.gamete_id, "OOCYTE_TEST_01");
        assert_eq!(fetched_gamete.tax_id, 9606);
        assert_eq!(fetched_gamete.mitochondrial_count, 220_000);
        assert_eq!(fetched_gamete.follicle_stage, "PreOvulatoryMII");

        // Longevity targets seeded during init_schema
        let targets = store.list_longevity_targets(9606).unwrap();
        assert_eq!(targets.len(), 10);
        let tert = store.get_longevity_target(9606, "TERT").unwrap().unwrap();
        assert_eq!(tert.category, "telomere");
        assert_eq!(tert.priority, 1);
        let sirt6 = store.get_longevity_target_by_symbol("SIRT6").unwrap().unwrap();
        assert_eq!(sirt6.category, "dna_repair");
        assert_eq!(sirt6.priority, 1);
    }
}
