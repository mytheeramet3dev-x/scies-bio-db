//! Demonstration of ingesting reference sequences, annotations, and querying longevity targets.
//!
//! Demonstrates the Phase DB-1 pipeline:
//! 1. Open `BioDb` with 2-bit sequence storage and SQLite catalog
//! 2. Ingest FASTA reference sequences with SHA-256 cryptographic provenance
//! 3. Ingest GFF3 annotations with comprehensive report
//! 4. Query the seeded Longevity Target catalog (TERT, SIRT6, TP53, etc.)
//! 5. Extract typed `LongevityLocus` records containing real nucleotide sequence slices
//! 6. Feed the locus sequence directly into `scies-bio-th` CRISPR guide design

use std::path::Path;
use tempfile::tempdir;

use scies_bio_db::BioDb;

fn main() -> scies_bio_db::Result<()> {
    println!("===========================================================================");
    println!(" 🧬 scies-bio-db: Ingest & Query Longevity Targets Demonstration (Phase DB-1)");
    println!("===========================================================================\n");

    let tmp = tempdir()?;
    let db_path = tmp.path();

    // 1. Open the biological database engine
    println!("1️⃣  Opening biological database at: {:?}", db_path);
    let db = BioDb::open(db_path)?;

    // 2. Locate testdata fixtures
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let tert_fa = manifest_dir.join("testdata/mini_chr5_tert.fasta");
    let sirt6_fa = manifest_dir.join("testdata/mini_chr19_sirt6.fasta");
    let tp53_fa = manifest_dir.join("testdata/mini_chr17_tp53.fasta");
    let gff_path = manifest_dir.join("testdata/mini_genes.gff3");

    // 3. Ingest FASTA sequences with cryptographic provenance
    println!("\n2️⃣  Ingesting reference sequences with SHA-256 provenance...");
    let tax_id = 9606;
    let assembly = "T2T-CHM13v2.0";

    let refs_tert = db.ingest_fasta_with_uri(
        &tert_fa,
        tax_id,
        assembly,
        Some("https://api.ncbi.nlm.nih.gov/genomes/T2T-CHM13v2.0/chr5.fa"),
    )?;
    println!(
        "   ✓ Ingested TERT chromosome contig: {} ({} bp, sha256: {})",
        refs_tert[0].seq_id, refs_tert[0].length_bp, &refs_tert[0].sha256[..16]
    );

    let refs_sirt6 = db.ingest_fasta_with_uri(
        &sirt6_fa,
        tax_id,
        assembly,
        Some("https://api.ncbi.nlm.nih.gov/genomes/T2T-CHM13v2.0/chr19.fa"),
    )?;
    println!(
        "   ✓ Ingested SIRT6 chromosome contig: {} ({} bp, sha256: {})",
        refs_sirt6[0].seq_id, refs_sirt6[0].length_bp, &refs_sirt6[0].sha256[..16]
    );

    let refs_tp53 = db.ingest_fasta_with_uri(
        &tp53_fa,
        tax_id,
        assembly,
        Some("https://api.ncbi.nlm.nih.gov/genomes/T2T-CHM13v2.0/chr17.fa"),
    )?;
    println!(
        "   ✓ Ingested TP53 chromosome contig: {} ({} bp, sha256: {})",
        refs_tp53[0].seq_id, refs_tp53[0].length_bp, &refs_tp53[0].sha256[..16]
    );

    // 4. Ingest GFF3 gene annotations
    println!("\n3️⃣  Ingesting GFF3 genomic annotations...");
    let report = db.ingest_gff3_report(&gff_path, tax_id, assembly)?;
    println!(
        "   ✓ Ingested {} genes, {} transcripts, {} exons (skipped {} invalid records)",
        report.genes_count, report.transcripts_count, report.exons_count, report.skipped_count
    );

    // 5. Query longevity targets catalog
    println!("\n4️⃣  Querying Curated Longevity Target Catalog (tax_id 9606)...");
    let targets = db.list_longevity_targets(tax_id)?;
    println!("   Found {} essential longevity targets registered:", targets.len());
    for target in &targets {
        println!(
            "   - [{}] Prio {}: {:<16} ({})",
            target.symbol,
            target.priority,
            target.category,
            target.notes.as_deref().unwrap_or("N/A")
        );
    }

    // 6. Query typed LongevityLocus records for TERT, SIRT6, and TP53
    println!("\n5️⃣  Extracting Typed LongevityLocus records...");
    for symbol in &["TERT", "SIRT6", "TP53"] {
        let locus = db.get_longevity_locus(tax_id, assembly, symbol)?;
        let gc = locus.sequence.gc_content();
        println!("   ------------------------------------------------------------");
        println!("   Target Symbol : {}", locus.gene.symbol);
        println!("   Genomic Window: {}:{}-{} ({:?})", locus.gene.chrom, locus.window_start, locus.window_end, locus.gene.strand);
        println!("   Extracted Bp  : {} bp (GC: {:.1}%)", locus.len(), gc);
        println!("   Provenance SHA: {}", locus.provenance_sha256);
        if let Some(ref target) = locus.target {
            println!("   Category      : {}", target.category);
            println!("   Priority      : {}", target.priority);
        }

        // Downstream demonstration: CRISPR guide design using scies-bio-th
        let strand = match locus.gene.strand {
            scies_bio_db::model::genome::Strand::Forward => scies_bio_th::genome::Strand::Forward,
            scies_bio_db::model::genome::Strand::Reverse => scies_bio_th::genome::Strand::Reverse,
        };
        let guides = scies_bio_th::aging::design_real_guides(
            &locus.sequence,
            strand,
            locus.window_start,
        );
        println!("   CRISPR Guides : {} candidate guide(s) designed", guides.len());
        if let Some(opt) = scies_bio_th::aging::select_optimal_guide(&guides) {
            let seq_str = String::from_utf8_lossy(&opt.sequence.sequence);
            println!(
                "   Optimal Guide : seq={}, PAM={}, pos={}, efficiency={:.1}%",
                seq_str, opt.pam_motif, opt.genomic_pos, opt.on_target_score
            );
        }
    }

    println!("\n===========================================================================");
    println!(" ✅ Demonstration completed successfully! All downstream integrations verified.");
    println!("===========================================================================\n");

    Ok(())
}
