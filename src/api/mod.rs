//! Top-level public API for `scies-bio-db`.
//!
//! [`BioDb`] is the single entry point that consumers of the library interact
//! with. It owns the storage backends and hands out typed query objects.

use std::path::Path;
use std::sync::Arc;

use crate::ingest::{
    BulkLoadReport, BulkLoader, FastaIngestor, GenomePackage, Gff3IngestReport, Gff3Ingestor,
    VcfIngestor,
};
use crate::model::gene::GeneFeature;
use crate::model::genome::{Assembly, SequenceRef};
use crate::model::longevity::{LongevityLocus, LongevityTarget};
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

    /// Clears any cached memory-mapped `.seq` file handles.
    pub fn clear_cache(&self) {
        self.seq_store.clear_cache();
    }

    /// Fetches a decoded sequence slice directly under `(tax_id, assembly, chr)` over half-open interval `[start, end)`.
    ///
    /// Preserves all IUPAC degenerate symbols and Uracil losslessly via `.seq v2`.
    pub fn fetch_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>> {
        self.seq_store.fetch(tax_id, assembly, chr, start, end)
    }

    /// Fetches a decoded DNA sequence slice over the 0-based half-open interval `[start, end)`
    /// and returns it as a [`scies_bio_th::Dna`] structure (Section 5 DB4).
    pub fn get_sequence(
        &self,
        tax_id: u32,
        assembly: &str,
        seq_id: &str,
        start: u64,
        end: u64,
    ) -> Result<scies_bio_th::Dna> {
        let raw = self.seq_store.fetch(tax_id, assembly, seq_id, start, end)?;
        scies_bio_th::Dna::new(&raw).map_err(crate::BioDbError::Bio)
    }

    /// Retrieves sequence metadata and cryptographic provenance for `(tax_id, assembly, seq_id)`.
    pub fn get_sequence_ref(
        &self,
        tax_id: u32,
        assembly: &str,
        seq_id: &str,
    ) -> Result<Option<SequenceRef>> {
        self.sql_store.get_sequence_ref(tax_id, assembly, seq_id)
    }

    /// Lists all sequences registered under `(tax_id, assembly)`.
    pub fn list_sequences(&self, tax_id: u32, assembly: &str) -> Result<Vec<SequenceRef>> {
        self.sql_store.list_sequences(tax_id, assembly)
    }

    /// Registers a reference assembly record in SQLite catalog.
    pub fn register_assembly(&self, assembly: &Assembly) -> Result<()> {
        self.sql_store.insert_assembly(assembly)
    }

    /// Retrieves reference assembly metadata.
    pub fn get_assembly(&self, tax_id: u32, assembly: &str) -> Result<Option<Assembly>> {
        self.sql_store.get_assembly(tax_id, assembly)
    }

    /// Lists all reference assemblies registered for `tax_id`.
    pub fn list_assemblies(&self, tax_id: u32) -> Result<Vec<Assembly>> {
        self.sql_store.list_assemblies(tax_id)
    }

    /// Ingests a FASTA file with an explicit `source_uri`, capturing cryptographic SHA-256 provenance.
    pub fn ingest_fasta_with_uri(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
        source_uri: Option<&str>,
    ) -> Result<Vec<SequenceRef>> {
        let ingestor = FastaIngestor::new(Arc::clone(&self.seq_store));
        ingestor.ingest_file_with_provenance(path, tax_id, assembly, source_uri, Some(&self.sql_store))
    }

    /// Verifies the cryptographic integrity of an ingested sequence by comparing its stored SHA-256
    /// hash in SQLite against the hash computed from reading the 2-bit sequence payload from disk.
    pub fn verify_sequence_integrity(
        &self,
        tax_id: u32,
        assembly: &str,
        seq_id: &str,
    ) -> Result<bool> {
        let seq_ref = self
            .sql_store
            .get_sequence_ref(tax_id, assembly, seq_id)?
            .ok_or_else(|| {
                crate::BioDbError::NotFound(format!(
                    "Sequence metadata not found for {tax_id}:{assembly}:{seq_id}"
                ))
            })?;

        let raw = self.seq_store.fetch(tax_id, assembly, seq_id, 0, seq_ref.length_bp)?;
        let computed = scies_bio_th::seq::compute_sha256_hex(&raw);

        if computed.eq_ignore_ascii_case(&seq_ref.sha256) {
            Ok(true)
        } else {
            Err(crate::BioDbError::CorruptData(format!(
                "Integrity failure for {tax_id}:{assembly}:{seq_id}: expected sha256 {}, computed {}",
                seq_ref.sha256, computed
            )))
        }
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

    /// Ingests a GFF3/GTF annotation file and returns a detailed [`Gff3IngestReport`] (Section 5 DB3).
    pub fn ingest_gff3_report(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
    ) -> Result<Gff3IngestReport> {
        let ingestor = Gff3Ingestor::new();
        ingestor.ingest_file_report(path, tax_id, assembly, &self.sql_store)
    }

    /// Finds a gene feature by symbol (or gene identifier) scoped to `(tax_id, assembly)` (Section 5 DB4).
    pub fn get_gene_by_symbol(
        &self,
        tax_id: u32,
        assembly: &str,
        symbol: &str,
    ) -> Result<GeneFeature> {
        self.sql_store
            .get_gene_feature_by_symbol(tax_id, assembly, symbol)?
            .ok_or_else(|| {
                crate::BioDbError::NotFound(format!(
                    "Gene symbol '{symbol}' not found in assembly '{assembly}' for tax_id {tax_id}"
                ))
            })
    }

    /// Lists all gene features overlapping 0-based half-open interval `[start, end)`
    /// on chromosome `seq_id` scoped to `(tax_id, assembly)` (Section 5 DB4).
    pub fn list_genes_in_region(
        &self,
        tax_id: u32,
        assembly: &str,
        seq_id: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<GeneFeature>> {
        if start > end {
            return Err(crate::BioDbError::InvalidCoordinate {
                chr: seq_id.to_string(),
                start,
                end,
            });
        }
        self.sql_store
            .list_gene_features_in_region(tax_id, assembly, seq_id, start, end)
    }

    // -----------------------------------------------------------------------
    // Longevity Targets & Locus Retrieval (Sprint DB-S3)
    // -----------------------------------------------------------------------

    /// Lists all registered longevity targets for `tax_id` sorted by priority.
    pub fn list_longevity_targets(&self, tax_id: u32) -> Result<Vec<LongevityTarget>> {
        self.sql_store.list_longevity_targets(tax_id)
    }

    /// Retrieves a longevity target by gene symbol (defaults to tax_id 9606 if multiple match).
    pub fn get_longevity_target(&self, symbol: &str) -> Result<Option<LongevityTarget>> {
        self.sql_store.get_longevity_target_by_symbol(symbol)
    }

    /// Retrieves a longevity target for a specific organism `tax_id` and `symbol`.
    pub fn get_longevity_target_for_organism(
        &self,
        tax_id: u32,
        symbol: &str,
    ) -> Result<Option<LongevityTarget>> {
        self.sql_store.get_longevity_target(tax_id, symbol)
    }

    /// Registers or updates a longevity target record in the catalog.
    pub fn register_longevity_target(&self, target: &LongevityTarget) -> Result<()> {
        self.sql_store.insert_longevity_target(target)
    }

    /// Fetches a fully resolved genomic locus for a longevity target gene in `(tax_id, assembly)`.
    ///
    /// Looks up the gene coordinate annotation, extracts the exact DNA sequence from the 2-bit
    /// sequence store over `[gene.start_bp, gene.end_bp)`, and attaches cryptographic SHA-256 provenance.
    pub fn get_longevity_locus(
        &self,
        tax_id: u32,
        assembly: &str,
        symbol: &str,
    ) -> Result<LongevityLocus> {
        self.get_longevity_locus_with_flanks(tax_id, assembly, symbol, 0, 0)
    }

    /// Fetches a complete genomic locus with customizable upstream and downstream coordinate flanks (in bp).
    ///
    /// Flanks expand the retrieved coordinate window:
    /// `[gene.start_bp.saturating_sub(upstream_bp), min(gene.end_bp + downstream_bp, seq_length)]`.
    pub fn get_longevity_locus_with_flanks(
        &self,
        tax_id: u32,
        assembly: &str,
        symbol: &str,
        upstream_bp: u64,
        downstream_bp: u64,
    ) -> Result<LongevityLocus> {
        let gene = self.get_gene_by_symbol(tax_id, assembly, symbol)?;

        let seq_ref = self
            .sql_store
            .get_sequence_ref(tax_id, assembly, &gene.chrom)?;

        let total_len = if let Some(ref sr) = seq_ref {
            sr.length_bp
        } else {
            self.seq_store.sequence_length(tax_id, assembly, &gene.chrom)?
        };

        let (fetch_start, fetch_end, window_start, window_end) = if total_len >= gene.end_bp {
            // Full chromosome or segment covering the absolute coordinates
            let ws = gene.start_bp.saturating_sub(upstream_bp);
            let we = (gene.end_bp + downstream_bp).min(total_len);
            (ws, we, ws, we)
        } else if total_len > 0 && gene.start_bp >= total_len {
            // Extracted mini contig whose origin (0) corresponds to gene.start_bp
            let ws = gene.start_bp;
            let we = gene.start_bp + total_len;
            (0, total_len, ws, we)
        } else {
            return Err(crate::BioDbError::InvalidCoordinate {
                chr: gene.chrom.clone(),
                start: gene.start_bp,
                end: gene.end_bp,
            });
        };

        if fetch_start >= fetch_end {
            return Err(crate::BioDbError::InvalidCoordinate {
                chr: gene.chrom.clone(),
                start: fetch_start,
                end: fetch_end,
            });
        }

        let raw = self
            .seq_store
            .fetch(tax_id, assembly, &gene.chrom, fetch_start, fetch_end)?;
        let sequence = scies_bio_th::Dna::new(&raw).map_err(crate::BioDbError::Bio)?;

        let target = self
            .sql_store
            .get_longevity_target(tax_id, symbol)?
            .or_else(|| self.sql_store.get_longevity_target_by_symbol(symbol).ok().flatten());

        let provenance_sha256 = seq_ref.map(|r| r.sha256).unwrap_or_default();

        Ok(LongevityLocus {
            gene,
            target,
            sequence,
            window_start,
            window_end,
            provenance_sha256,
            tax_id,
            assembly: assembly.to_string(),
        })
    }

    /// Ingests a VCF variant file into the metadata store under `(tax_id, assembly)`.
    ///
    /// Returns the number of variants ingested.
    pub fn ingest_vcf(&self, path: &Path, tax_id: u32, assembly: &str) -> Result<u64> {
        let ingestor = VcfIngestor::new();
        ingestor.ingest_file(path, tax_id, assembly, &self.sql_store)
    }

    /// Bulk-loads an entire NCBI or Ensembl [`GenomePackage`] into the database atomically.
    ///
    /// Ingests sequences into `.seq v2` format (leveraging multi-core parallelism when the `parallel` feature is active),
    /// imports GFF3 features in transactions, records VCF variants, and returns a comprehensive [`BulkLoadReport`].
    pub fn bulk_load_package(&self, package: &GenomePackage) -> Result<BulkLoadReport> {
        let loader = BulkLoader::new(Arc::clone(&self.seq_store), Arc::clone(&self.sql_store));
        loader.load(package)
    }

    /// Automatically scans a stored chromosome sequence for Gardiner-Garden & Frommer CpG islands
    /// using `scies-bio-th::epigenome` and registers them as [`crate::model::EpigeneticFeature`] records in `BioDb`.
    pub fn index_chromosome_cpg_islands(
        &self,
        tax_id: u32,
        assembly: &str,
        chr: &str,
    ) -> Result<usize> {
        let dna = self.genome().fetch_full_sequence(tax_id, assembly, chr)?;
        let islands =
            scies_bio_th::epigenome::methylation::find_cpg_islands(&dna, 200, 50, 0.50, 0.60);
        let count = islands.len();
        for (idx, isl) in islands.into_iter().enumerate() {
            let feat = crate::model::EpigeneticFeature {
                feature_id: format!("CGI:{}:{}:{}", assembly, chr, idx + 1),
                tax_id,
                assembly_name: assembly.to_string(),
                chr: chr.to_string(),
                start: isl.start as u64,
                end: isl.end as u64,
                track_type: crate::model::EpigeneticTrackType::CpgIsland,
                mark_label: format!("CpG_Island_OE_{:.2}", isl.obs_exp_ratio),
                signal_value: isl.obs_exp_ratio,
            };
            self.sql_store.insert_epigenetic_feature(&feat)?;
        }
        Ok(count)
    }

    /// Loads all stored [`crate::model::RegulatoryEdgeRecord`] interactions for `tax_id`
    /// and constructs an in-memory [`scies_bio_th::regulatory::GeneRegulatoryNetwork`].
    pub fn load_grn(
        &self,
        tax_id: u32,
    ) -> Result<scies_bio_th::regulatory::GeneRegulatoryNetwork> {
        let edges = self.sql_store.list_regulatory_edges(tax_id)?;
        let mut grn = scies_bio_th::regulatory::GeneRegulatoryNetwork::new();
        for e in edges {
            grn.add_regulation(&e.tf_symbol, &e.target_symbol, e.interaction, e.weight);
        }
        Ok(grn)
    }

    /// Queries a stored [`crate::model::MultiOmicsTrackRecord`] for `(tax_id, assembly, sample_id, gene_symbol)`,
    /// computes its proximal enhancer ABC score via `scies_bio_th::regulatory`, and fuses all layers
    /// via `scies_bio_th::multiomics::fuse_multiomics_locus`.
    pub fn load_and_fuse_locus(
        &self,
        tax_id: u32,
        assembly: &str,
        sample_id: &str,
        gene_symbol: &str,
    ) -> Result<scies_bio_th::multiomics::MultiOmicsFusionResult> {
        let track = self
            .sql_store
            .get_multiomics_track(tax_id, assembly, sample_id, gene_symbol)?
            .ok_or_else(|| {
                crate::error::BioDbError::NotFound(format!(
                    "Multi-omics track not found for {}:{}:{}",
                    assembly, sample_id, gene_symbol
                ))
            })?;

        let enh = scies_bio_th::regulatory::score_enhancer_activity(
            track.h3k27ac_chip,
            track.h3k27ac_input,
            track.atac_accessibility,
            track.enhancer_distance_bp as usize,
        );

        let fusion_input = scies_bio_th::multiomics::OmicsLocusInput {
            gene_symbol: track.gene_symbol,
            rna_tpm: track.rna_tpm,
            atac_accessibility: track.atac_accessibility,
            promoter_methylation_beta: track.promoter_methylation_beta,
            enhancer_abc_score: enh.abc_score,
        };

        Ok(scies_bio_th::multiomics::fuse_multiomics_locus(&fusion_input))
    }

    /// Loads spatial transcriptomics cells and decoded mRNA spots for `dataset_id` from SQLite
    /// and assigns each transcript to its nearest cell within `max_radius_um` via `scies_bio_th::imaging`.
    pub fn assign_spatial_transcripts(
        &self,
        dataset_id: &str,
        max_radius_um: f64,
    ) -> Result<Vec<scies_bio_th::imaging::SpatialCellAssignment>> {
        let cell_recs = self.sql_store.list_spatial_cells(dataset_id)?;
        let spot_recs = self.sql_store.list_spatial_transcripts(dataset_id)?;

        let cells: Vec<scies_bio_th::imaging::CellCentroid> = cell_recs
            .into_iter()
            .map(|c| scies_bio_th::imaging::CellCentroid {
                cell_id: c.cell_id,
                x: c.x_um,
                y: c.y_um,
            })
            .collect();

        let spots: Vec<scies_bio_th::imaging::TranscriptSpot> = spot_recs
            .into_iter()
            .map(|s| scies_bio_th::imaging::TranscriptSpot {
                gene_symbol: s.gene_symbol,
                x: s.x_um,
                y: s.y_um,
            })
            .collect();

        Ok(scies_bio_th::imaging::assign_transcripts_to_cells(
            &cells,
            &spots,
            max_radius_um,
        ))
    }

    /// Loads single-cell states for `dataset_id` from SQLite and computes a Diffusion Pseudotime
    /// trajectory rooted at `root_cell_id` via `scies_bio_th::singlecell`.
    pub fn compute_atlas_pseudotime(
        &self,
        dataset_id: &str,
        root_cell_id: &str,
    ) -> Result<scies_bio_th::singlecell::PseudotimeTrajectory> {
        let recs = self.sql_store.list_single_cell_records(dataset_id)?;
        let nodes: Vec<scies_bio_th::singlecell::CellStateNode> = recs
            .into_iter()
            .map(|r| scies_bio_th::singlecell::CellStateNode {
                cell_id: r.cell_id,
                expression: r.expression,
            })
            .collect();
        Ok(scies_bio_th::singlecell::compute_diffusion_pseudotime(
            &nodes,
            root_cell_id,
        )?)
    }

    /// Queries a stored [`crate::model::PharmaCompoundRecord`] by `compound_id` and simulates
    /// a repeated dosing regimen via `scies_bio_th::pharma::simulate_pkpd_regimen`.
    pub fn simulate_compound_regimen(
        &self,
        compound_id: &str,
        dose_um: f64,
        interval_hours: usize,
        num_doses: usize,
    ) -> Result<scies_bio_th::pharma::PkPdRegimenReport> {
        let cmpd = self
            .sql_store
            .get_pharma_compound(compound_id)?
            .ok_or_else(|| {
                crate::error::BioDbError::NotFound(format!(
                    "Pharmacokinetic compound '{}' not found",
                    compound_id
                ))
            })?;
        let params = scies_bio_th::pharma::PkPdParameters {
            compound_name: cmpd.compound_name,
            ec50_um: cmpd.ec50_um,
            td50_um: cmpd.td50_um,
            hill_coefficient: cmpd.hill_coefficient,
            half_life_hours: cmpd.half_life_hours,
            bioavailability: cmpd.bioavailability,
        };
        Ok(scies_bio_th::pharma::simulate_pkpd_regimen(
            &params,
            dose_um,
            interval_hours,
            num_doses,
        ))
    }

    /// Synthesizes an MII oocyte using the 3-stage `rOvary` organoid folliculogenesis model
    /// (`scies_bio_th::ivg::organoid::simulate_rovary_folliculogenesis`) from a donor genome and stores the
    /// resulting [`crate::model::SyntheticGameteRecord`] in SQLite.
    pub fn synthesize_and_store_oocyte(
        &self,
        gamete_id: &str,
        tax_id: u32,
        donor_dna: &scies_bio_th::Dna,
        config: &scies_bio_th::ivg::organoid::OvarianOrganoidConfig,
        meiotic_seed: u64,
    ) -> Result<crate::model::SyntheticGameteRecord> {
        let pgclc = scies_bio_th::ivg::PgclcCulture::induce(
            donor_dna,
            tax_id,
            scies_bio_th::ivg::InductionProtocol::Saitou2018,
            8,
        )?;
        let rep = scies_bio_th::ivg::organoid::simulate_rovary_folliculogenesis(
            &pgclc,
            config,
            meiotic_seed,
        )?;
        let record = crate::model::SyntheticGameteRecord {
            gamete_id: gamete_id.to_string(),
            tax_id,
            gamete_type: "Oocyte".to_string(),
            follicle_stage: format!("{:?}", rep.final_stage),
            diameter_um: rep.oocyte_diameter_um,
            mitochondrial_count: rep.oocyte.quality.mitochondrial_count,
            imprint_fidelity: rep.oocyte.quality.imprint_fidelity,
            aneuploidy_risk: rep.oocyte.quality.aneuploidy_risk,
            overall_score: rep.oocyte.quality.overall_score,
        };
        self.sql_store.insert_synthetic_gamete(&record)?;
        Ok(record)
    }

    /// Retrieves a [`crate::model::SyntheticGameteRecord`] by `gamete_id`.
    pub fn get_synthetic_gamete(
        &self,
        gamete_id: &str,
    ) -> Result<Option<crate::model::SyntheticGameteRecord>> {
        self.sql_store.get_synthetic_gamete(gamete_id)
    }

    /// Executes the high-fidelity longevity genome design pipeline (`scies_bio_th::aging::design_high_fidelity_longevity_package`)
    /// using stored reference sequences in `BioDb`, then scans every proposed sgRNA across ALL chromosomes
    /// in `assembly` via [`GenomeQuery::scan_assembly_crispr_offtargets`] and updates the report's off-target
    /// scores and oncogenic risk assessment with assembly-wide cross-chromosome specificity.
    pub fn design_and_audit_assembly_longevity_package(
        &self,
        tax_id: u32,
        assembly: &str,
        subject: &scies_bio_th::aging::HighFidelityAgeParams,
        strategy: scies_bio_th::aging::LongevityStrategy,
        gene_chr_map: &[(&str, &str)],
        max_mismatches: usize,
    ) -> Result<scies_bio_th::aging::HighFidelityDesignReport> {
        let mut owned_loci: Vec<(String, scies_bio_th::Dna)> =
            Vec::with_capacity(gene_chr_map.len());
        for &(gene_sym, chr_name) in gene_chr_map {
            let dna = self.genome().fetch_full_sequence(tax_id, assembly, chr_name)?;
            owned_loci.push((gene_sym.to_string(), dna));
        }

        let loci_refs: Vec<(&str, &scies_bio_th::Dna)> = owned_loci
            .iter()
            .map(|(sym, dna)| (sym.as_str(), dna))
            .collect();

        let mut report = scies_bio_th::aging::design_high_fidelity_longevity_package(
            subject,
            strategy,
            Some(&loci_refs),
        );

        // Re-evaluate each guide's off-target score across all stored chromosomes in the assembly
        for edit in &mut report.proposed_edits {
            if let Some(ref mut guide) = edit.recommended_sgrna {
                let scan = self.genome().scan_assembly_crispr_offtargets(
                    tax_id,
                    assembly,
                    &guide.sequence,
                    max_mismatches,
                )?;
                guide.offtarget_cfd_score = scan.specificity_score;
                edit.offtarget_score = (100.0 - scan.specificity_score).clamp(1.0, 50.0);
            }
        }

        let telomerase_effective = match strategy {
            scies_bio_th::aging::LongevityStrategy::Conservative => 0.0,
            scies_bio_th::aging::LongevityStrategy::Balanced => 0.80,
            scies_bio_th::aging::LongevityStrategy::Aggressive => 0.95,
        };
        let offtarget_fractions: Vec<f64> = report
            .proposed_edits
            .iter()
            .map(|e| e.offtarget_score / 100.0)
            .collect();
        report.oncogenic_risk = scies_bio_th::aging::evaluate_oncogenic_risk(
            telomerase_effective,
            &offtarget_fractions,
            report.metrics.projected_shelterin_capping,
            &report.oncogenic_risk.safeguards,
        );

        Ok(report)
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

        // 7. Verify direct BioDb::fetch_sequence and IUPAC/RNA raw fetch
        let direct_slice = db.fetch_sequence(9606, "GRCh38", "chr1", 0, 4).unwrap();
        assert_eq!(&direct_slice, b"ATGC");

        let iupac_seq = b"ACGTUNRYSWKMBDHV";
        db.store_sequence(9606, "GRCh38", "chr_iupac", iupac_seq)
            .unwrap();
        let raw_fetched = db
            .genome()
            .fetch_raw_sequence(9606, "GRCh38", "chr_iupac", 0, iupac_seq.len() as u64)
            .unwrap();
        assert_eq!(&raw_fetched, iupac_seq);
    }

    #[test]
    fn test_bulk_package_and_indices() {
        let dir = tempdir().unwrap();
        let db = BioDb::open(dir.path()).unwrap();

        let ecoli = Organism {
            tax_id: 562,
            scientific_name: "Escherichia coli".to_string(),
            common_name: Some("E. coli".to_string()),
            lineage: Lineage {
                domain: "Bacteria".into(),
                kingdom: "Bacteria".into(),
                phylum: "Pseudomonadota".into(),
                class: "Gammaproteobacteria".into(),
                order: "Enterobacterales".into(),
                family: "Enterobacteriaceae".into(),
                genus: "Escherichia".into(),
                species: "coli".into(),
            },
        };
        let genome = ReferenceGenome {
            assembly_name: "ASM584v2".to_string(),
            tax_id: 562,
            status: crate::model::genome::AssemblyStatus::Primary,
            release_year: 2014,
        };

        // Create sample fasta
        let mut fa = NamedTempFile::new().unwrap();
        writeln!(fa, ">chr1 E. coli genome").unwrap();
        writeln!(
            fa,
            "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC"
        )
        .unwrap();
        fa.flush().unwrap();

        // Create sample gff3
        let mut gff = NamedTempFile::new().unwrap();
        writeln!(gff, "##gff-version 3").unwrap();
        writeln!(
            gff,
            "chr1\tRefSeq\tgene\t1\t30\t.\t+\t.\tID=gene:lacZ;Name=lacZ;biotype=protein_coding"
        )
        .unwrap();
        writeln!(
            gff,
            "chr1\tRefSeq\tmRNA\t1\t30\t.\t+\t.\tID=transcript:lacZ_t1;Parent=gene:lacZ;biotype=mRNA"
        )
        .unwrap();
        writeln!(
            gff,
            "chr1\tRefSeq\texon\t1\t30\t.\t+\t.\tID=exon:lacZ_e1;Parent=transcript:lacZ_t1;exon_number=1"
        )
        .unwrap();
        gff.flush().unwrap();

        let package = GenomePackage::new(ecoli, genome)
            .add_fasta(fa.path())
            .add_gff3(gff.path());

        let report = db.bulk_load_package(&package).unwrap();
        assert_eq!(report.chromosomes_ingested, 1);
        assert!(report.total_bases_ingested > 50);
        assert_eq!(report.features_ingested, 3);

        // Verify Interval Index
        let interval_idx = db
            .genes()
            .build_chromosome_interval_index(562, "ASM584v2", "chr1")
            .unwrap();
        assert_eq!(interval_idx.len(), 1);
        let hits = interval_idx.find_overlapping("chr1", 5, 15);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].payload_id, "lacZ");

        // Verify K-mer similarity
        let sim = db
            .genome()
            .estimate_similarity(562, "ASM584v2", "chr1", b"ATGACCATGATTAC", 6)
            .unwrap();
        assert!(sim > 0.0);

        // Verify CRISPR off-target scanner on stored chromosome
        let guide = scies_bio_th::Rna::new(b"AUGACCAUGAUUACGCCAAG").unwrap();
        let offtargets = db
            .genome()
            .scan_crispr_offtargets(562, "ASM584v2", "chr1", &guide, 5)
            .unwrap();
        // Should return a valid list without error
        assert!(offtargets.len() <= 10);

        // Verify Epigenomic CpG island auto-indexing, ProteinStructure, and DatasetManifest
        let cpg_rich = b"CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
        db.store_sequence(562, "ASM584v2", "chr_cpg", cpg_rich).unwrap();
        let cpg_count = db
            .index_chromosome_cpg_islands(562, "ASM584v2", "chr_cpg")
            .unwrap();
        assert!(cpg_count >= 1);
        let epi_features = db
            .sql_store()
            .find_epigenetic_features_by_region(562, "ASM584v2", "chr_cpg", 0, 250)
            .unwrap();
        assert_eq!(epi_features.len(), cpg_count);

        let struct_rec = crate::model::epigenome::ProteinStructureRecord {
            structure_id: "AF-P00722-F1".to_string(),
            gene_symbol: "lacZ".to_string(),
            tax_id: 562,
            resolution_angstrom: 1.85,
            residue_count: 1024,
            pocket_volume_a3: 412.5,
        };
        db.sql_store().insert_protein_structure(&struct_rec).unwrap();
        let fetched_structs = db.sql_store().find_protein_structures_by_gene(562, "lacZ").unwrap();
        assert_eq!(fetched_structs.len(), 1);
        assert!((fetched_structs[0].pocket_volume_a3 - 412.5).abs() < 1e-6);

        let manifest = crate::model::epigenome::DatasetManifest {
            manifest_id: "refseq-ecoli-2026".to_string(),
            tax_id: 562,
            assembly_name: "ASM584v2".to_string(),
            provider: "NCBI RefSeq".to_string(),
            source_sha256: "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855".to_string(),
            ingested_at: "2026-11-03T10:00:00Z".to_string(),
        };
        db.sql_store().insert_dataset_manifest(&manifest).unwrap();
        let fetched_manifest = db.sql_store().get_dataset_manifest("refseq-ecoli-2026").unwrap().unwrap();
        assert_eq!(fetched_manifest.provider, "NCBI RefSeq");

        // Verify v0.4 Regulatory GRN, Multi-Omics Fusion, Promoter Scanner & Spatial Omics
        db.sql_store()
            .insert_regulatory_edge(&crate::model::RegulatoryEdgeRecord {
                edge_id: "ECOLI_CRP_LACZ".to_string(),
                tax_id: 562,
                tf_symbol: "CRP".to_string(),
                target_symbol: "lacZ".to_string(),
                interaction: scies_bio_th::regulatory::RegulationType::Activation,
                weight: 0.90,
                evidence_source: "RegulonDB".to_string(),
            })
            .unwrap();
        let grn = db.load_grn(562).unwrap();
        assert_eq!(grn.topological_order().unwrap(), vec!["CRP", "lacZ"]);

        db.sql_store()
            .insert_multiomics_track(&crate::model::MultiOmicsTrackRecord {
                tax_id: 562,
                assembly_name: "ASM584v2".to_string(),
                sample_id: "WT_LogPhase".to_string(),
                gene_symbol: "lacZ".to_string(),
                rna_tpm: 512.0,
                atac_accessibility: 0.88,
                promoter_methylation_beta: 0.08,
                h3k27ac_chip: 25.0,
                h3k27ac_input: 1.5,
                enhancer_distance_bp: 3000,
            })
            .unwrap();
        let fused = db
            .load_and_fuse_locus(562, "ASM584v2", "WT_LogPhase", "lacZ")
            .unwrap();
        assert!(fused.unified_activity_score > 0.75);

        let prom_seq = b"GGGCGCCTATAAAAGGCCCCGGGCGGGGTCAGTTT";
        db.store_sequence(562, "ASM584v2", "chr_prom", prom_seq)
            .unwrap();
        let prom_scan = db
            .genome()
            .scan_promoter_elements(562, "ASM584v2", "chr_prom", 0, prom_seq.len() as u64)
            .unwrap();
        assert!(!prom_scan.core_elements.is_empty());
        assert!(!prom_scan.tfbs_hits.is_empty());

        db.sql_store()
            .insert_spatial_cell(&crate::model::SpatialCellRecord {
                dataset_id: "MERFISH_01".to_string(),
                cell_id: "Cell_A".to_string(),
                x_um: 10.0,
                y_um: 10.0,
            })
            .unwrap();
        db.sql_store()
            .insert_spatial_transcript(&crate::model::SpatialTranscriptRecord {
                dataset_id: "MERFISH_01".to_string(),
                spot_id: "Spot_1".to_string(),
                gene_symbol: "lacZ".to_string(),
                x_um: 11.0,
                y_um: 9.5,
            })
            .unwrap();
        let spatial_assigned = db.assign_spatial_transcripts("MERFISH_01", 5.0).unwrap();
        assert_eq!(spatial_assigned.len(), 1);
        assert_eq!(spatial_assigned[0].total_transcripts, 1);

        // Verify v0.5 Single-Cell Atlas Pseudotime & Pharma PK/PD Regimen Simulation
        db.sql_store()
            .insert_single_cell_record(&crate::model::SingleCellAtlasRecord {
                dataset_id: "ATLAS_REPROG".to_string(),
                cell_id: "Cell_Root".to_string(),
                cell_type: "Somatic".to_string(),
                expression: vec![0.0, 0.0, 0.0],
            })
            .unwrap();
        db.sql_store()
            .insert_single_cell_record(&crate::model::SingleCellAtlasRecord {
                dataset_id: "ATLAS_REPROG".to_string(),
                cell_id: "Cell_Mid".to_string(),
                cell_type: "PartialReprog".to_string(),
                expression: vec![1.0, 1.0, 0.0],
            })
            .unwrap();
        let traj = db
            .compute_atlas_pseudotime("ATLAS_REPROG", "Cell_Root")
            .unwrap();
        assert_eq!(traj.points.len(), 2);
        assert_eq!(traj.points[0].cell_id, "Cell_Root");

        db.sql_store()
            .insert_pharma_compound(&crate::model::PharmaCompoundRecord {
                compound_id: "CMPD_FISETIN".to_string(),
                compound_name: "Fisetin".to_string(),
                target_gene: "BCL2L1".to_string(),
                ec50_um: 2.0,
                td50_um: 25.0,
                hill_coefficient: 1.8,
                half_life_hours: 8.0,
                bioavailability: 0.85,
            })
            .unwrap();
        let pkpd_rep = db
            .simulate_compound_regimen("CMPD_FISETIN", 5.0, 12, 4)
            .unwrap();
        assert!(pkpd_rep.within_therapeutic_window);

        // Verify v1.0.0 Synthetic Gamete Generation & Storage
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

        let organoid_cfg = scies_bio_th::ivg::organoid::OvarianOrganoidConfig::standard_rovary();
        let donor_dna = scies_bio_th::Dna::new(b"ATGCGATCGATCGATCGATCGATCGATCGA").unwrap();
        let stored_gamete = db
            .synthesize_and_store_oocyte("OOCYTE_HUMAN_EXP_01", 9606, &donor_dna, &organoid_cfg, 42)
            .unwrap();
        assert_eq!(stored_gamete.gamete_id, "OOCYTE_HUMAN_EXP_01");
        assert_eq!(stored_gamete.gamete_type, "Oocyte");
        assert_eq!(stored_gamete.follicle_stage, "PreOvulatoryMII");
        assert!(stored_gamete.diameter_um >= 105.0);
        assert!(stored_gamete.mitochondrial_count >= 100_000);
        assert!(stored_gamete.overall_score >= 80.0);

        let fetched_gamete = db
            .get_synthetic_gamete("OOCYTE_HUMAN_EXP_01")
            .unwrap()
            .unwrap();
        assert_eq!(fetched_gamete.gamete_id, "OOCYTE_HUMAN_EXP_01");
        assert_eq!(fetched_gamete.tax_id, 9606);
        assert_eq!(fetched_gamete.follicle_stage, "PreOvulatoryMII");
    }

    #[test]
    fn test_multi_chromosome_grch38_assembly_crispr_offtarget_scan() {
        let dir = tempdir().unwrap();
        let db = BioDb::open(dir.path()).unwrap();

        let human = Organism {
            tax_id: 9606,
            scientific_name: "Homo sapiens".to_string(),
            common_name: Some("Human".to_string()),
            lineage: Lineage {
                domain: "Eukaryota".into(),
                kingdom: "Animalia".into(),
                phylum: "Chordata".into(),
                class: "Mammalia".into(),
                order: "Primates".into(),
                family: "Hominidae".into(),
                genus: "Homo".into(),
                species: "sapiens".into(),
            },
        };
        db.register_organism(&human).unwrap();
        db.register_reference_genome(&ReferenceGenome {
            assembly_name: "GRCh38.p14".to_string(),
            tax_id: 9606,
            status: crate::model::genome::AssemblyStatus::Primary,
            release_year: 2022,
        })
        .unwrap();

        // Ingest real GRCh38 loci across 5 chromosomes from scies-bio-th/testdata
        let testdata_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("scies-bio-th")
            .join("testdata");
        for file in &[
            "t2t_chr17_tp53.fasta",
            "t2t_chr19_sirt6.fasta",
            "t2t_chr9_cdkn2a.fasta",
            "t2t_chr6_foxo3.fasta",
            "t2t_chr11_atm.fasta",
        ] {
            db.ingest_fasta(&testdata_dir.join(file), 9606, "GRCh38.p14")
                .unwrap();
        }

        // Design a real guide on chr17_tp53_locus and scan across all 5 chromosomes in GRCh38.p14
        let tp53_dna = db
            .genome()
            .fetch_full_sequence(9606, "GRCh38.p14", "chr17_tp53_locus")
            .unwrap();
        let guides = scies_bio_th::aging::design_real_guides(
            &tp53_dna,
            scies_bio_th::genome::Strand::Reverse,
            7_687_001,
        );
        let optimal = scies_bio_th::aging::select_optimal_guide(&guides).unwrap();

        let report = db
            .genome()
            .scan_assembly_crispr_offtargets(9606, "GRCh38.p14", &optimal.sequence, 4)
            .unwrap();

        assert_eq!(report.chromosomes_scanned, 5);
        assert_eq!(report.total_bp_scanned, 4_400);
        assert_eq!(report.exact_matches, 1);
        assert_eq!(report.hits[0].chr, "chr17_tp53_locus");
        assert_eq!(report.hits[0].mismatches, 0);
        assert!((report.hits[0].cfd_score - 100.0).abs() < 1e-6);
        assert!(report.specificity_score > 50.0 && report.specificity_score <= 100.0);

        // Verify full assembly-wide longevity design + cross-chromosome off-target audit
        db.ingest_fasta(&testdata_dir.join("t2t_chr1_parp1.fasta"), 9606, "GRCh38.p14")
            .unwrap();
        db.ingest_fasta(&testdata_dir.join("t2t_chr5_telomere_tert.fasta"), 9606, "GRCh38.p14")
            .unwrap();

        let subject = scies_bio_th::aging::HighFidelityAgeParams {
            chronological_age: 58.0,
            cpg_methylation_fraction: 0.52,
            telomere_length_bp: 6_100.0,
            gamma_h2ax_foci: 1.9,
            sasp_burden: 0.38,
        };
        let gene_chr_map = [
            ("TP53", "chr17_tp53_locus"),
            ("CDKN2A", "chr9_cdkn2a_locus"),
            ("TERT", "chr5_tert_promoter"),
            ("SIRT6", "chr19_sirt6_locus"),
            ("PARP1", "chr1_parp1_locus"),
            ("FOXO3", "chr6_foxo3_locus"),
        ];
        let full_report = db
            .design_and_audit_assembly_longevity_package(
                9606,
                "GRCh38.p14",
                &subject,
                scies_bio_th::aging::LongevityStrategy::Balanced,
                &gene_chr_map,
                4,
            )
            .unwrap();

        assert!(full_report
            .proposed_edits
            .iter()
            .all(|e| e.sequence_provenance == scies_bio_th::aging::SequenceProvenance::RealReference));
        assert!(full_report.oncogenic_risk.is_acceptable_safety_profile);
    }

    #[test]
    fn test_sprint_db_s1_storage_fasta_and_provenance() {
        let dir = tempdir().unwrap();
        let db = BioDb::open(dir.path()).unwrap();

        // 1. Ingest mini FASTA fixtures
        let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let tert_fa = manifest_dir.join("testdata/mini_chr5_tert.fasta");
        let sirt6_fa = manifest_dir.join("testdata/mini_chr19_sirt6.fasta");
        let tp53_fa = manifest_dir.join("testdata/mini_chr17_tp53.fasta");

        let chrs5 = db.ingest_fasta(&tert_fa, 9606, "T2T-CHM13v2.0").unwrap();
        assert_eq!(chrs5, vec!["chr5"]);

        let chrs19 = db.ingest_fasta(&sirt6_fa, 9606, "T2T-CHM13v2.0").unwrap();
        assert_eq!(chrs19, vec!["chr19"]);

        let chrs17 = db.ingest_fasta(&tp53_fa, 9606, "T2T-CHM13v2.0").unwrap();
        assert_eq!(chrs17, vec!["chr17"]);

        // 2. Query get_sequence(tax_id, assembly, seq_id, start, end) -> Result<Dna>
        let dna5 = db.get_sequence(9606, "T2T-CHM13v2.0", "chr5", 0, 20).unwrap();
        assert_eq!(dna5.sequence, b"CCGTGCGTCCCAGGGCACGC");
        assert_eq!(dna5.sequence.len(), 20);

        // Substring inside sequence
        let sub19 = db.get_sequence(9606, "T2T-CHM13v2.0", "chr19", 0, 15).unwrap();
        assert_eq!(sub19.sequence, b"GGAGGGGGCAGCCGA");

        let sub17 = db.get_sequence(9606, "T2T-CHM13v2.0", "chr17", 0, 12).unwrap();
        assert_eq!(sub17.sequence, b"TTCGAGATCCAA");

        // 3. Coordinate validation
        let invalid_coord = db.get_sequence(9606, "T2T-CHM13v2.0", "chr5", 100, 50);
        assert!(matches!(invalid_coord, Err(crate::BioDbError::InvalidCoordinate { .. })));

        let out_of_bounds = db.get_sequence(9606, "T2T-CHM13v2.0", "chr5", 0, 999_999);
        assert!(matches!(out_of_bounds, Err(crate::BioDbError::InvalidCoordinate { .. })));

        // 4. Provenance tracking via sequences & assemblies catalog
        let sref5 = db.get_sequence_ref(9606, "T2T-CHM13v2.0", "chr5").unwrap().expect("chr5 ref exists");
        assert_eq!(sref5.seq_id, "chr5");
        assert_eq!(sref5.length_bp, 800);
        assert!(!sref5.sha256.is_empty());
        assert_eq!(sref5.path, "seq/9606/T2T-CHM13v2.0/chr5.seq");

        let seq_list = db.list_sequences(9606, "T2T-CHM13v2.0").unwrap();
        assert_eq!(seq_list.len(), 3);

        let assemblies = db.list_assemblies(9606).unwrap();
        assert_eq!(assemblies.len(), 1);
        assert_eq!(assemblies[0].assembly_name, "T2T-CHM13v2.0");

        // 5. Sequence integrity check
        let intact = db.verify_sequence_integrity(9606, "T2T-CHM13v2.0", "chr5").unwrap();
        assert!(intact);

        // 6. Detect corruption
        db.clear_cache();
        let seq_file_path = dir.path().join("seq/9606/T2T-CHM13v2.0/chr5.seq");
        let mut bytes = std::fs::read(&seq_file_path).unwrap();
        let corrupt_idx = bytes.len() - 1;
        bytes[corrupt_idx] ^= 0x55;
        std::fs::write(&seq_file_path, bytes).unwrap();

        let verify_result = db.verify_sequence_integrity(9606, "T2T-CHM13v2.0", "chr5");
        assert!(verify_result.is_err(), "Corrupted sequence file must fail verification");

        // 7. Re-ingest updated sequence -> hash updates and integrity passes
        db.clear_cache();
        let new_fa_path = dir.path().join("chr5_v2.fasta");
        std::fs::write(&new_fa_path, b">chr5\nACGTACGTACGT\n").unwrap();
        db.ingest_fasta(&new_fa_path, 9606, "T2T-CHM13v2.0").unwrap();
        let sref_v2 = db.get_sequence_ref(9606, "T2T-CHM13v2.0", "chr5").unwrap().unwrap();
        assert_ne!(sref_v2.sha256, sref5.sha256);
        assert_eq!(sref_v2.length_bp, 12);
        assert!(db.verify_sequence_integrity(9606, "T2T-CHM13v2.0", "chr5").unwrap());
    }

    #[test]
    fn test_sprint_db_s2_annotation_and_gene_queries() {
        let dir = tempdir().unwrap();
        let db = BioDb::open(dir.path()).unwrap();

        let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let gff_path = manifest_dir.join("testdata/mini_genes.gff3");

        // 1. Ingest mini GFF3 with report
        let report = db
            .ingest_gff3_report(&gff_path, 9606, "T2T-CHM13v2.0")
            .unwrap();
        assert_eq!(report.genes_count, 3);
        assert_eq!(report.transcripts_count, 3);
        assert_eq!(report.exons_count, 3);
        assert_eq!(report.skipped_count, 1); // 1 malformed record skipped
        assert_eq!(report.total_inserted, 9);

        // 2. get_gene_by_symbol tests (TERT, SIRT6, TP53)
        let tert = db
            .get_gene_by_symbol(9606, "T2T-CHM13v2.0", "TERT")
            .unwrap();
        assert_eq!(tert.symbol, "TERT");
        assert_eq!(tert.chrom, "chr5");
        assert_eq!(tert.start_bp, 1294800);
        assert_eq!(tert.end_bp, 1295600);
        assert_eq!(tert.strand, crate::model::genome::Strand::Forward);
        assert_eq!(tert.feature_type, "gene");

        let sirt6 = db
            .get_gene_by_symbol(9606, "T2T-CHM13v2.0", "SIRT6")
            .unwrap();
        assert_eq!(sirt6.symbol, "SIRT6");
        assert_eq!(sirt6.chrom, "chr19");
        assert_eq!(sirt6.start_bp, 4182000);
        assert_eq!(sirt6.end_bp, 4182800);
        assert_eq!(sirt6.strand, crate::model::genome::Strand::Forward);

        let tp53 = db
            .get_gene_by_symbol(9606, "T2T-CHM13v2.0", "TP53")
            .unwrap();
        assert_eq!(tp53.symbol, "TP53");
        assert_eq!(tp53.chrom, "chr17");
        assert_eq!(tp53.start_bp, 7687000);
        assert_eq!(tp53.end_bp, 7688000);
        assert_eq!(tp53.strand, crate::model::genome::Strand::Reverse);

        // 3. Error case for missing symbol
        let missing = db.get_gene_by_symbol(9606, "T2T-CHM13v2.0", "UNKNOWN_GENE");
        assert!(matches!(missing, Err(crate::BioDbError::NotFound(_))));

        // 4. Region query list_genes_in_region
        let region_matches = db
            .list_genes_in_region(9606, "T2T-CHM13v2.0", "chr5", 1294000, 1296000)
            .unwrap();
        assert_eq!(region_matches.len(), 1);
        assert_eq!(region_matches[0].symbol, "TERT");

        let empty_region = db
            .list_genes_in_region(9606, "T2T-CHM13v2.0", "chr5", 0, 1000)
            .unwrap();
        assert!(empty_region.is_empty());

        let invalid_coords = db.list_genes_in_region(9606, "T2T-CHM13v2.0", "chr5", 5000, 4000);
        assert!(matches!(
            invalid_coords,
            Err(crate::BioDbError::InvalidCoordinate { .. })
        ));
    }

    #[test]
    fn test_sprint_db_s3_longevity_catalog_and_locus() {
        let dir = tempdir().unwrap();
        let db = BioDb::open(dir.path()).unwrap();

        // 1. Seeded longevity targets check (tax_id 9606)
        let targets = db.list_longevity_targets(9606).unwrap();
        assert_eq!(targets.len(), 10, "10 essential longevity genes must be seeded");

        // Verify priorities: priority 1 comes before priority 2
        let p1: Vec<&str> = targets
            .iter()
            .filter(|t| t.priority == 1)
            .map(|t| t.symbol.as_str())
            .collect();
        assert_eq!(p1, vec!["SIRT6", "TERC", "TERT", "TP53"]);

        let p2: Vec<&str> = targets
            .iter()
            .filter(|t| t.priority == 2)
            .map(|t| t.symbol.as_str())
            .collect();
        assert_eq!(p2, vec!["ATM", "CDKN2A", "FOXO3", "KLOTHO", "LMNA", "PARP1"]);

        // 2. Query target by symbol
        let tert_target = db.get_longevity_target("TERT").unwrap().expect("TERT target exists");
        assert_eq!(tert_target.symbol, "TERT");
        assert_eq!(tert_target.category, "telomere");
        assert_eq!(tert_target.priority, 1);
        assert_eq!(tert_target.default_assembly, Some("T2T-CHM13v2.0".to_string()));

        let sirt6_target = db.get_longevity_target("SIRT6").unwrap().expect("SIRT6 target exists");
        assert_eq!(sirt6_target.category, "dna_repair");

        let tp53_target = db.get_longevity_target_for_organism(9606, "TP53").unwrap().expect("TP53 target exists");
        assert_eq!(tp53_target.category, "tumor_suppressor");

        // 3. Register custom target
        let custom_target = LongevityTarget {
            symbol: "FOXM1".to_string(),
            tax_id: 9606,
            priority: 3,
            category: "stemness".to_string(),
            default_assembly: Some("T2T-CHM13v2.0".to_string()),
            notes: Some("Forkhead box M1; pluripotency and mitotic regulator".to_string()),
        };
        db.register_longevity_target(&custom_target).unwrap();
        let fetched_custom = db.get_longevity_target("FOXM1").unwrap().unwrap();
        assert_eq!(fetched_custom.category, "stemness");
        assert_eq!(fetched_custom.priority, 3);

        // 4. Ingest fixtures
        let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let tert_fa = manifest_dir.join("testdata/mini_chr5_tert.fasta");
        let sirt6_fa = manifest_dir.join("testdata/mini_chr19_sirt6.fasta");
        let tp53_fa = manifest_dir.join("testdata/mini_chr17_tp53.fasta");
        let gff_path = manifest_dir.join("testdata/mini_genes.gff3");

        db.ingest_fasta(&tert_fa, 9606, "T2T-CHM13v2.0").unwrap();
        db.ingest_fasta(&sirt6_fa, 9606, "T2T-CHM13v2.0").unwrap();
        db.ingest_fasta(&tp53_fa, 9606, "T2T-CHM13v2.0").unwrap();
        db.ingest_gff3_report(&gff_path, 9606, "T2T-CHM13v2.0").unwrap();

        // 5. Fetch full LongevityLocus for TERT
        let locus = db.get_longevity_locus(9606, "T2T-CHM13v2.0", "TERT").unwrap();
        assert_eq!(locus.gene.symbol, "TERT");
        assert_eq!(locus.gene.chrom, "chr5");
        assert_eq!(locus.gene.start_bp, 1294800);
        assert_eq!(locus.gene.end_bp, 1295600);
        assert_eq!(locus.window_start, 1294800);
        assert_eq!(locus.window_end, 1295600);
        assert_eq!(locus.len(), 800);
        assert!(!locus.is_empty());
        assert!(!locus.provenance_sha256.is_empty(), "Provenance SHA-256 must be populated");
        assert_eq!(locus.tax_id, 9606);
        assert_eq!(locus.assembly, "T2T-CHM13v2.0");

        let target_meta = locus.target.as_ref().expect("Target metadata must be attached");
        assert_eq!(target_meta.symbol, "TERT");
        assert_eq!(target_meta.category, "telomere");

        // Verify sequence starts with expected prefix
        assert!(locus.sequence.sequence.starts_with(b"CCGTGCGTCC"));

        // 6. Fetch SIRT6 and TP53 loci
        let sirt6_locus = db.get_longevity_locus(9606, "T2T-CHM13v2.0", "SIRT6").unwrap();
        assert_eq!(sirt6_locus.gene.symbol, "SIRT6");
        assert_eq!(sirt6_locus.len(), 800);
        assert!(sirt6_locus.sequence.sequence.starts_with(b"GGAGGGGGCA"));

        let tp53_locus = db.get_longevity_locus(9606, "T2T-CHM13v2.0", "TP53").unwrap();
        assert_eq!(tp53_locus.gene.symbol, "TP53");
        assert_eq!(tp53_locus.len(), 1000);
        assert!(tp53_locus.sequence.sequence.starts_with(b"TTCGAGATCC"));

        // 7. Error cases
        let missing_gene = db.get_longevity_locus(9606, "T2T-CHM13v2.0", "NONEXISTENT_GENE");
        assert!(matches!(missing_gene, Err(crate::BioDbError::NotFound(_))));
    }
}
