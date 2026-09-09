//! Top-level public API for `scies-bio-db`.
//!
//! [`BioDb`] is the single entry point that consumers of the library interact
//! with. It owns the storage backends and hands out typed query objects.

use std::path::Path;
use std::sync::Arc;

use crate::ingest::{FastaIngestor, Gff3Ingestor, VcfIngestor};
use crate::model::taxonomy::Organism;
use crate::query::{GeneQuery, GenomeQuery, ProteinQuery, VariantQuery};
use crate::storage::{SeqStore, SqliteStore};
use crate::Result;

/// The main entry point for `scies-bio-db`.
///
/// Open a database directory with [`BioDb::open`] and use the typed accessor
/// methods to query biological data.
///
/// # Example
///
/// ```no_run
/// use scies_bio_db::BioDb;
/// use std::path::Path;
///
/// let db = BioDb::open(Path::new("./bio_db")).unwrap();
/// let dna = db.genome().fetch_sequence(9606, "GRCh38", "chr1", 0, 100).unwrap();
/// ```
pub struct BioDb {
    /// Shared 2-bit sequence storage.
    seq_store: Arc<SeqStore>,
    /// Shared SQLite metadata store.
    sql_store: Arc<SqliteStore>,
}

impl BioDb {
    /// Opens or creates a [`BioDb`] rooted at `db_dir`.
    ///
    /// `db_dir` is created if it does not exist. The sequence store lives in
    /// `<db_dir>/seq/` and the SQLite file at `<db_dir>/meta.db`.
    pub fn open(db_dir: &Path) -> Result<Self> {
        std::fs::create_dir_all(db_dir)?;
        let seq_dir = db_dir.join("seq");
        let meta_path = db_dir.join("meta.db");

        let seq_store = SeqStore::new(&seq_dir)?;
        let sql_store = SqliteStore::open(&meta_path)?;
        sql_store.init_schema()?;

        Ok(Self {
            seq_store: Arc::new(seq_store),
            sql_store: Arc::new(sql_store),
        })
    }

    /// Returns a reference to the underlying sequence store.
    pub fn seq_store(&self) -> &Arc<SeqStore> {
        &self.seq_store
    }

    /// Returns a reference to the underlying SQLite store.
    pub fn sql_store(&self) -> &Arc<SqliteStore> {
        &self.sql_store
    }

    /// Stores a raw DNA sequence directly under `(tax_id, assembly, chr)`.
    pub fn store_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        sequence: &[u8],
    ) -> Result<()> {
        self.seq_store.store(tax_id, assembly, chr, sequence)
    }

    /// Registers or updates an organism in the database.
    pub fn register_organism(&self, organism: &Organism) -> Result<()> {
        self.sql_store.insert_organism(organism)
    }

    /// Registers a reference genome assembly in the database.
    pub fn register_reference_genome(
        &self,
        genome: &crate::model::genome::ReferenceGenome,
    ) -> Result<()> {
        self.sql_store.insert_reference_genome(genome)
    }

    /// Retrieves an organism by taxonomy ID.
    pub fn get_organism(&self, tax_id: u32) -> Result<Option<Organism>> {
        self.sql_store.get_organism(tax_id)
    }

    /// Lists all registered organisms in the database.
    pub fn list_organisms(&self) -> Result<Vec<Organism>> {
        self.sql_store.list_organisms()
    }

    /// Returns a [`GenomeQuery`] handle for chromosome and sequence queries.
    pub fn genome(&self) -> GenomeQuery {
        GenomeQuery::new(Arc::clone(&self.seq_store), Arc::clone(&self.sql_store))
    }

    /// Returns a [`GeneQuery`] handle for gene and transcript queries.
    pub fn genes(&self) -> GeneQuery {
        GeneQuery::new(Arc::clone(&self.sql_store))
    }

    /// Returns a [`ProteinQuery`] handle for protein entry and domain queries.
    pub fn proteins(&self) -> ProteinQuery {
        ProteinQuery::new(Arc::clone(&self.sql_store))
    }

    /// Returns a [`VariantQuery`] handle for variant record queries.
    pub fn variants(&self) -> VariantQuery {
        VariantQuery::new(Arc::clone(&self.sql_store))
    }

    /// Ingests a FASTA file, encoding sequences into the sequence store under `(tax_id, assembly)`.
    ///
    /// Returns the list of chromosome names that were ingested.
    pub fn ingest_fasta(&self, path: &Path, tax_id: u32, assembly: &str) -> Result<Vec<String>> {
        let ingestor = FastaIngestor::new(Arc::clone(&self.seq_store));
        ingestor.ingest_file(path, tax_id, assembly, Some(&self.sql_store))
    }

    /// Ingests a GFF3 annotation file into the metadata store under `(tax_id, assembly)`.
    ///
    /// Returns the number of features ingested.
    pub fn ingest_gff3(&self, path: &Path, tax_id: u32, assembly: &str) -> Result<u64> {
        let ingestor = Gff3Ingestor::new();
        ingestor.ingest_file(path, tax_id, assembly, &self.sql_store)
    }

    /// Ingests a VCF variant file into the metadata store under `(tax_id, assembly)`.
    ///
    /// Returns the number of variants ingested.
    pub fn ingest_vcf(&self, path: &Path, tax_id: u32, assembly: &str) -> Result<u64> {
        let ingestor = VcfIngestor::new();
        ingestor.ingest_file(path, tax_id, assembly, &self.sql_store)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::genome::ReferenceGenome;
    use crate::model::taxonomy::Lineage;
    use std::io::Write;
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_biodb_end_to_end() {
        let dir = tempdir().unwrap();
        let db = BioDb::open(dir.path()).unwrap();

        // 1. Register Human (TaxID: 9606)
        let human = Organism {
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
        db.register_organism(&human).unwrap();

        // Register Mouse (TaxID: 10090)
        let mouse = Organism {
            tax_id: 10090,
            scientific_name: "Mus musculus".to_string(),
            common_name: Some("Mouse".to_string()),
            lineage: Lineage {
                domain: "Eukaryota".to_string(),
                kingdom: "Animalia".to_string(),
                phylum: "Chordata".to_string(),
                class: "Mammalia".to_string(),
                order: "Rodentia".to_string(),
                family: "Muridae".to_string(),
                genus: "Mus".to_string(),
                species: "musculus".to_string(),
            },
        };
        db.register_organism(&mouse).unwrap();

        let orgs = db.list_organisms().unwrap();
        assert_eq!(orgs.len(), 2);

        // 2. Register Reference Genomes
        db.register_reference_genome(&ReferenceGenome {
            assembly_name: "GRCh38".to_string(),
            tax_id: 9606,
            status: crate::model::genome::AssemblyStatus::Primary,
            release_year: 2013,
        })
        .unwrap();

        db.register_reference_genome(&ReferenceGenome {
            assembly_name: "GRCm39".to_string(),
            tax_id: 10090,
            status: crate::model::genome::AssemblyStatus::Primary,
            release_year: 2020,
        })
        .unwrap();

        // 3. Ingest FASTA sequences for both organisms
        let mut human_fa = NamedTempFile::new().unwrap();
        writeln!(human_fa, ">chr1 Human Chromosome 1").unwrap();
        writeln!(human_fa, "ATGCGATCGATC").unwrap();
        human_fa.flush().unwrap();

        let mut mouse_fa = NamedTempFile::new().unwrap();
        writeln!(mouse_fa, ">chr1 Mouse Chromosome 1").unwrap();
        writeln!(mouse_fa, "GGCCGGCCGGCC").unwrap();
        mouse_fa.flush().unwrap();

        db.ingest_fasta(human_fa.path(), 9606, "GRCh38").unwrap();
        db.ingest_fasta(mouse_fa.path(), 10090, "GRCm39").unwrap();

        // 4. Fetch and verify sequences are isolated by tax_id / assembly
        let human_dna = db
            .genome()
            .fetch_sequence(9606, "GRCh38", "chr1", 0, 12)
            .unwrap();
        assert_eq!(human_dna.sequence, b"ATGCGATCGATC");

        let mouse_dna = db
            .genome()
            .fetch_sequence(10090, "GRCm39", "chr1", 0, 12)
            .unwrap();
        assert_eq!(mouse_dna.sequence, b"GGCCGGCCGGCC");

        // 5. Verify chromosome listing from SQLite
        let human_chrs = db.genome().list_chromosomes("GRCh38").unwrap();
        assert_eq!(human_chrs.len(), 1);
        assert_eq!(human_chrs[0].name, "chr1");
        assert_eq!(human_chrs[0].length, 12);

        // 6. Test Multi-Species & Multi-Assembly Scoped Gene Queries
        let human_gene = crate::model::gene::Gene {
            gene_id: "ENSG00000141510".to_string(),
            assembly_name: "GRCh38".to_string(),
            symbol: "TP53".to_string(),
            name: "tumor protein p53".to_string(),
            biotype: crate::model::gene::Biotype::ProteinCoding,
            chr: "chr17".to_string(),
            start: 7668401,
            end: 7687550,
            strand: crate::model::genome::Strand::Reverse,
            tax_id: 9606,
        };
        db.sql_store().insert_gene(&human_gene).unwrap();

        let mouse_gene = crate::model::gene::Gene {
            gene_id: "ENSMUSG00000059552".to_string(),
            assembly_name: "GRCm39".to_string(),
            symbol: "Trp53".to_string(),
            name: "transformation related protein 53".to_string(),
            biotype: crate::model::gene::Biotype::ProteinCoding,
            chr: "chr11".to_string(),
            start: 69482126,
            end: 69494399,
            strand: crate::model::genome::Strand::Reverse,
            tax_id: 10090,
        };
        db.sql_store().insert_gene(&mouse_gene).unwrap();

        // Cross-species isolation test with identical symbol
        let human_shared = crate::model::gene::Gene {
            gene_id: "ENSG_SHARED".to_string(),
            assembly_name: "GRCh38".to_string(),
            symbol: "SHARED_GENE".to_string(),
            name: "shared gene human".to_string(),
            biotype: crate::model::gene::Biotype::ProteinCoding,
            chr: "chr1".to_string(),
            start: 100,
            end: 500,
            strand: crate::model::genome::Strand::Forward,
            tax_id: 9606,
        };
        let mouse_shared = crate::model::gene::Gene {
            gene_id: "ENSMUS_SHARED".to_string(),
            assembly_name: "GRCm39".to_string(),
            symbol: "SHARED_GENE".to_string(),
            name: "shared gene mouse".to_string(),
            biotype: crate::model::gene::Biotype::ProteinCoding,
            chr: "chr1".to_string(),
            start: 100,
            end: 500,
            strand: crate::model::genome::Strand::Forward,
            tax_id: 10090,
        };
        db.sql_store().insert_gene(&human_shared).unwrap();
        db.sql_store().insert_gene(&mouse_shared).unwrap();

        let fetched_h = db
            .genes()
            .find_by_symbol_scoped(9606, "SHARED_GENE")
            .unwrap();
        assert_eq!(fetched_h.gene_id, "ENSG_SHARED");
        assert_eq!(fetched_h.assembly_name, "GRCh38");

        let fetched_m = db
            .genes()
            .find_by_symbol_scoped(10090, "SHARED_GENE")
            .unwrap();
        assert_eq!(fetched_m.gene_id, "ENSMUS_SHARED");
        assert_eq!(fetched_m.assembly_name, "GRCm39");

        // Region isolation check
        let h_region = db
            .genes()
            .find_by_region(9606, "GRCh38", "chr1", 50, 600)
            .unwrap();
        assert_eq!(h_region.len(), 1);
        assert_eq!(h_region[0].gene_id, "ENSG_SHARED");

        let m_region = db
            .genes()
            .find_by_region(10090, "GRCm39", "chr1", 50, 600)
            .unwrap();
        assert_eq!(m_region.len(), 1);
        assert_eq!(m_region[0].gene_id, "ENSMUS_SHARED");
    }
}
