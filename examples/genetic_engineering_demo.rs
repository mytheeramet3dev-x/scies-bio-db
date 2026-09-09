//! Genetic Engineering & CRISPR Genome Editing Demonstration
//! Using `scies-bio-db` and `scies-bio-th`.

use scies_bio_db::BioDb;
use scies_bio_th::edit::crispr::{find_cas12a_pam, find_spcas9_pam, generate_sgrna_candidates};
use scies_bio_th::edit::offtarget::calculate_cfd_score;
use scies_bio_th::edit::primer::{check_gc_clamp, detect_dimer, detect_hairpin, Primer};
use scies_bio_th::seq::Dna;
use std::path::Path;

fn main() -> scies_bio_db::Result<()> {
    println!("=========================================================================");
    println!(" ✂️  scies-bio-th & scies-bio-db: Genetic Engineering & CRISPR Pipeline");
    println!("=========================================================================\n");

    let db_path = Path::new("./demo_biodb");
    if !db_path.exists() {
        eprintln!("Please run real_dna_demo first to initialize the database!");
        return Ok(());
    }

    let db = BioDb::open(db_path)?;

    // =======================================================================
    // 🧬 WORKFLOW 1: CRISPR-Cas9 Target Design on Human Beta-Globin (HBB)
    //    Application: Sickle Cell Anemia Gene Editing (Exon 1 / Codon 6)
    // =======================================================================
    println!("-------------------------------------------------------------------------");
    println!("🎯 WORKFLOW 1: CRISPR-Cas9 Guide RNA Design on Human Hemoglobin Beta (HBB)");
    println!("-------------------------------------------------------------------------");

    let hbb_dna = db
        .genome()
        .fetch_full_sequence(9606, "GRCh38.p14", "chr11_HBB")?;
    println!(
        "📖 Retrieved Human HBB gene ({} bp) from scies-bio-db",
        hbb_dna.sequence.len()
    );

    // 1. Find SpCas9 PAM sites (NGG)
    let spcas9_pams = find_spcas9_pam(&hbb_dna);
    println!("🔎 SpCas9 (NGG) PAM sites found: {}", spcas9_pams.len());
    for (i, pam) in spcas9_pams.iter().take(5).enumerate() {
        println!(
            "   PAM #{}: Position {} -> {}",
            i + 1,
            pam.index,
            String::from_utf8_lossy(&pam.sequence)
        );
    }

    // 2. Find Cas12a PAM sites (TTTV)
    let cas12a_pams = find_cas12a_pam(&hbb_dna);
    println!("\n🔎 Cas12a (TTTV) PAM sites found: {}", cas12a_pams.len());
    for (i, pam) in cas12a_pams.iter().take(3).enumerate() {
        println!(
            "   PAM #{}: Position {} -> {}",
            i + 1,
            pam.index,
            String::from_utf8_lossy(&pam.sequence)
        );
    }

    // 3. Generate sgRNA 20nt guide candidates for SpCas9
    let sgrna_candidates = generate_sgrna_candidates(&hbb_dna);
    println!(
        "\n🧬 sgRNA Candidates (20nt protospacer before NGG PAM): {}",
        sgrna_candidates.len()
    );
    for (i, cand) in sgrna_candidates.iter().take(4).enumerate() {
        println!(
            "   [Guide #{}] PAM at pos {:3} | 5'-{}-3'",
            i + 1,
            cand.pam_index,
            String::from_utf8_lossy(&cand.sequence.sequence)
        );
    }

    // 4. Evaluate Off-Target Risk with CFD Score
    if let Some(top_guide) = sgrna_candidates.first() {
        println!("\n🛡️  Evaluating Guide #1 Specificity & Off-Target Risk...");
        let guide_rna = &top_guide.sequence;
        println!(
            "   Guide RNA: 5'-{}-3'",
            String::from_utf8_lossy(&guide_rna.sequence)
        );

        // Exact match target (on-target)
        let on_target_dna =
            Dna::new(&hbb_dna.sequence[top_guide.pam_index - 20..top_guide.pam_index])
                .map_err(scies_bio_db::BioDbError::Bio)?;
        let on_target_score = calculate_cfd_score(guide_rna, &on_target_dna);
        println!(
            "   👉 On-target match score: {:.1}% (mismatches: {})",
            on_target_score.score, on_target_score.mismatches
        );

        // Hypothetical off-target with 2 mismatches outside seed
        let mut off1_bytes = on_target_dna.sequence.clone();
        off1_bytes[0] = b'T'; // 5' end mismatch
        off1_bytes[2] = b'A';
        let off1_dna = Dna::new(&off1_bytes).map_err(scies_bio_db::BioDbError::Bio)?;
        let off1_score = calculate_cfd_score(guide_rna, &off1_dna);
        println!(
            "   👉 Off-target (2 non-seed mismatches) CFD score: {:.2}% (mismatches: {})",
            off1_score.score, off1_score.mismatches
        );

        // Hypothetical off-target with mismatch in seed region (near PAM)
        let mut off2_bytes = on_target_dna.sequence.clone();
        off2_bytes[18] = b'G'; // Seed region mismatch
        let off2_dna = Dna::new(&off2_bytes).map_err(scies_bio_db::BioDbError::Bio)?;
        let off2_score = calculate_cfd_score(guide_rna, &off2_dna);
        println!("   👉 Off-target (1 seed-region mismatch) CFD score: {:.2}% (mismatches: {}) -> Heavily penalized!", off2_score.score, off2_score.mismatches);
    }

    // =======================================================================
    // 🧪 WORKFLOW 2: PCR Primer Design for Human Insulin (INS) Cloning
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🧪 WORKFLOW 2: PCR Primer Design to Clone Human Insulin into E. coli");
    println!("-------------------------------------------------------------------------");

    let ins_dna = db
        .genome()
        .fetch_full_sequence(9606, "GRCh38.p14", "chr11_INS")?;
    println!(
        "📖 Retrieved Human Insulin CDS ({} bp)",
        ins_dna.sequence.len()
    );

    // Design Forward Primer (first 20 bp)
    let fwd_seq = Dna::new(&ins_dna.sequence[..20]).map_err(scies_bio_db::BioDbError::Bio)?;
    let fwd_primer = Primer::new(fwd_seq);

    // Design Reverse Primer (reverse complement of last 20 bp)
    let rev_template = Dna::new(&ins_dna.sequence[ins_dna.sequence.len() - 20..])
        .map_err(scies_bio_db::BioDbError::Bio)?;
    let rev_primer_seq = rev_template.reverse_complement();
    let rev_primer = Primer::new(rev_primer_seq);

    println!(
        "   🔬 Forward Primer: 5'-{}-3'",
        String::from_utf8_lossy(&fwd_primer.sequence.sequence)
    );
    println!(
        "      - Length: {} bp | Tm: {:.1}°C | GC: {:.1}% | GC Clamp: {}",
        fwd_primer.sequence.sequence.len(),
        fwd_primer.tm,
        fwd_primer.gc_percent,
        check_gc_clamp(&fwd_primer.sequence)
    );

    println!(
        "   🔬 Reverse Primer: 5'-{}-3'",
        String::from_utf8_lossy(&rev_primer.sequence.sequence)
    );
    println!(
        "      - Length: {} bp | Tm: {:.1}°C | GC: {:.1}% | GC Clamp: {}",
        rev_primer.sequence.sequence.len(),
        rev_primer.tm,
        rev_primer.gc_percent,
        check_gc_clamp(&rev_primer.sequence)
    );

    let fwd_hairpin = detect_hairpin(&fwd_primer.sequence);
    let rev_hairpin = detect_hairpin(&rev_primer.sequence);
    let dimer_check = detect_dimer(&fwd_primer.sequence, &rev_primer.sequence);

    println!(
        "      - Hairpin Secondary Structure: Forward = {}, Reverse = {}",
        if fwd_hairpin {
            "⚠️ Hairpin"
        } else {
            "✅ Clean"
        },
        if rev_hairpin {
            "⚠️ Hairpin"
        } else {
            "✅ Clean"
        }
    );
    println!(
        "      - Primer Dimer Risk: {}",
        if dimer_check {
            "⚠️ WARNING: Dimer detected"
        } else {
            "✅ PASSED (No strong dimer)"
        }
    );

    // =======================================================================
    // 🧬 WORKFLOW 3: Recombinant Plasmid Construction (In Silico Cloning)
    //    Inserting Human Insulin into Bacterial Expression Vector (pET-28a style)
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🏭 WORKFLOW 3: In Silico Recombinant Expression Vector Construction");
    println!("-------------------------------------------------------------------------");

    // Bacterial expression cassette parts:
    let t7_promoter = b"TAATACGACTCACTATAGGG"; // T7 promoter
    let rbs_site = b"AAGGAGATATACAT"; // Ribosome Binding Site (Shine-Dalgarno)
    let his_tag = b"ATGCACCACCACCACCACCAC"; // 6xHis-tag for purification
    let ins_insert = &ins_dna.sequence; // Human INS gene
    let t7_term = b"CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG"; // T7 Terminator

    let mut recombinant_plasmid_cassette = Vec::new();
    recombinant_plasmid_cassette.extend_from_slice(t7_promoter);
    recombinant_plasmid_cassette.extend_from_slice(b"GAATTGTTATCCGCTCACAATTC"); // Lac operator
    recombinant_plasmid_cassette.extend_from_slice(rbs_site);
    recombinant_plasmid_cassette.extend_from_slice(his_tag);
    recombinant_plasmid_cassette.extend_from_slice(ins_insert);
    recombinant_plasmid_cassette.extend_from_slice(t7_term);

    let recombinant_dna =
        Dna::new(&recombinant_plasmid_cassette).map_err(scies_bio_db::BioDbError::Bio)?;

    println!("   ✨ Recombinant Cassette Constructed:");
    println!("      [T7 Promoter] -> [lac Operator] -> [RBS] -> [6xHis-Tag] -> [Human INS] -> [T7 Terminator]");
    println!(
        "      Total construct size: {} bp",
        recombinant_dna.sequence.len()
    );
    println!("      GC Content: {:.2}%", recombinant_dna.gc_content());

    // Save Recombinant Vector to scies-bio-db under E. coli engineered strain!
    println!("\n💾 Storing Engineered Construct into scies-bio-db (E. coli expression host)...");
    db.seq_store().store(
        562,
        "Recombinant_pET28_INS",
        "plasmid_insert",
        &recombinant_dna.sequence,
    )?;
    println!("   ✅ Successfully saved under: seq/562/Recombinant_pET28_INS/plasmid_insert.seq");

    // Fetch back and verify
    let stored_construct =
        db.genome()
            .fetch_full_sequence(562, "Recombinant_pET28_INS", "plasmid_insert")?;
    assert_eq!(stored_construct.sequence, recombinant_dna.sequence);
    println!(
        "   ✅ Verified stored construct integrity: 100% match ({} bp)",
        stored_construct.sequence.len()
    );

    println!("\n=========================================================================");
    println!(" 🎉 Genetic Engineering, CRISPR & Recombinant DNA Pipeline Complete!");
    println!("=========================================================================\n");

    Ok(())
}
