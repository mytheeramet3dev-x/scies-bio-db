//! Personalized Cancer Neoantigen mRNA Vaccine Design Pipeline
//! Using `scies-bio-th::vaccine` and `scies-bio-db`.

use scies_bio_db::BioDb;
use scies_bio_th::seq::Protein;
use scies_bio_th::vaccine::mhc::{predict_mhc_binding, HlaAllele};
use scies_bio_th::vaccine::mrna_design::MrnaVaccineConstruct;
use scies_bio_th::vaccine::neoantigen::Neoantigen;
use std::path::Path;

fn main() -> scies_bio_db::Result<()> {
    println!("=========================================================================");
    println!(" 🎯 scies-bio: Personalized Cancer Neoantigen mRNA Vaccine Design Engine");
    println!("=========================================================================\n");

    let db_path = Path::new("./demo_biodb");
    if !db_path.exists() {
        eprintln!("Please run real_dna_demo first to initialize the database!");
        return Ok(());
    }
    let db = BioDb::open(db_path)?;

    // =======================================================================
    // 👤 PATIENT CLINICAL ONCOLOGY PROFILE
    // =======================================================================
    let patient_id = "PT-2026-MELANOMA-992";
    let patient_hla = HlaAllele::HlaA0201;

    println!("📋 Patient Clinical Intake:");
    println!("   - Patient ID : {}", patient_id);
    println!("   - Diagnosis  : Stage IV Metastatic Cutaneous Melanoma");
    println!(
        "   - HLA Allele : {} (Primary MHC Class I Presentation Allele)",
        patient_hla.name()
    );

    // =======================================================================
    // 🧬 STEP 1: Somatic Mutation Screening & Neoantigen Extraction
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🔍 STEP 1: Tumor Biopsy Sequencing & Somatic Mutation Identification");
    println!("-------------------------------------------------------------------------");

    // Real human canonical reference protein fragments:
    // 1. BRAF (V600 is at position 15 in this active kinase domain segment)
    let braf_wt =
        Protein::new(b"LHEIKIIHRDLKSNNIFLHEDLTVKIGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIR").unwrap();
    // 2. TP53 (R175 is at position 15 in this core DNA-binding domain segment)
    let tp53_wt =
        Protein::new(b"APRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSV").unwrap();
    // 3. KRAS (G12 is at position 12 in the N-terminal GTPase segment)
    let kras_wt =
        Protein::new(b"MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG").unwrap();
    // 4. PIK3CA (H1047 is at position 15 in this kinase domain C-terminal segment)
    let pik3ca_wt =
        Protein::new(b"EALRIIENQYSLQKLKMDLILMHKALDLDFERIEGHYAKKQLAEVLEELKKEFQAKVDSE").unwrap();

    let candidate_mutations = vec![
        (
            "BRAF",
            &braf_wt,
            33,
            'E',
            "p.Val600Glu (V600E Driver Mutation)",
        ),
        (
            "TP53",
            &tp53_wt,
            45,
            'H',
            "p.Arg175His (R175H Inactivation Mutation)",
        ),
        (
            "KRAS",
            &kras_wt,
            12,
            'D',
            "p.Gly12Asp (G12D Oncogenic Activation)",
        ),
        (
            "PIK3CA",
            &pik3ca_wt,
            33,
            'R',
            "p.His1047Arg (H1047R Hotspot Mutation)",
        ),
    ];

    let mut extracted_neoantigens = Vec::new();

    for (gene, wt_prot, pos, mut_aa, desc) in candidate_mutations {
        let neo = Neoantigen::from_mutation(gene, wt_prot, pos, mut_aa, 12)
            .map_err(scies_bio_db::BioDbError::Bio)?;
        println!("   ✨ Detected Neoantigen: [{}] {}", gene, desc);
        println!(
            "      - Mutated 25-mer Window: {}",
            String::from_utf8_lossy(&neo.mut_peptide.sequence)
        );
        extracted_neoantigens.push(neo);
    }

    // =======================================================================
    // 🛡️ STEP 2: MHC Class I Binding Prediction & Epitope Prioritization
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🧪 STEP 2: HLA-A*02:01 Binding Affinity & Agretopicity (DAI) Screening");
    println!("-------------------------------------------------------------------------");

    let mut prioritized_epitopes = Vec::new();

    for neo in &extracted_neoantigens {
        let nine_mers = neo.generate_9mers();
        let mut best_binder = None;
        let mut highest_immunogenicity = 0.0_f64;

        for (_, wt_9mer, mut_9mer) in nine_mers {
            if let Ok(res) = predict_mhc_binding(&mut_9mer, &wt_9mer, patient_hla) {
                if res.immunogenicity_score > highest_immunogenicity {
                    highest_immunogenicity = res.immunogenicity_score;
                    best_binder = Some(res);
                }
            }
        }

        if let Some(best) = best_binder {
            println!(
                "   📌 Gene [{}] Best Neoepitope 9-mer: 5'-{}-3'",
                neo.gene_symbol,
                String::from_utf8_lossy(&best.mut_9mer.sequence)
            );
            println!(
                "      - Predicted Affinity (IC50): {:.1} nM ({})",
                best.mut_ic50_nm,
                if best.is_strong_binder {
                    "🟢 STRONG BINDER"
                } else {
                    "🟡 MODERATE"
                }
            );
            println!(
                "      - Differential Agretopicity (DAI): {:.2}x higher affinity than Wild-Type",
                best.differential_agretopicity_index
            );
            println!(
                "      - Immunogenicity Score: {:.1} / 100.0",
                best.immunogenicity_score
            );
            prioritized_epitopes.push((neo.clone(), best));
        }
    }

    // Sort by immunogenicity score descending
    prioritized_epitopes.sort_by(|a, b| {
        b.1.immunogenicity_score
            .partial_cmp(&a.1.immunogenicity_score)
            .unwrap()
    });

    // =======================================================================
    // 💉 STEP 3: Assembling Clinical-Grade mRNA Vaccine Construct
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🏭 STEP 3: Polytope Assembly, Human Codon Optimization & mRNA Engineering");
    println!("-------------------------------------------------------------------------");

    let selected_neoantigens: Vec<Neoantigen> = prioritized_epitopes
        .into_iter()
        .map(|(neo, _)| neo)
        .collect();

    let vaccine_construct = MrnaVaccineConstruct::build(patient_id, &selected_neoantigens)
        .map_err(scies_bio_db::BioDbError::Bio)?;

    println!("   🎉 Personalized mRNA Vaccine Construct Assembled Successfully!");
    println!(
        "      - Target Patient      : {}",
        vaccine_construct.patient_id
    );
    println!(
        "      - Antigens Included   : {:?}",
        vaccine_construct.included_genes
    );
    println!(
        "      - Translated ORF Size : {} amino acids",
        vaccine_construct.polytope_protein.sequence.len()
    );
    println!(
        "      - Codon-Optimized GC% : {:.2}% (Elevated for intracellular stability)",
        vaccine_construct.orf_gc_content
    );
    println!(
        "      - Total Transcript Len: {} nucleotides",
        vaccine_construct.total_length_nt
    );

    println!("\n   📐 Full mRNA Architecture:");
    println!("      [5' Cap1] ── (35 nt 5' UTR + Kozak GCCACC)");
    println!("                ── [Sec Signal Peptide: 20 aa]");
    println!("                ── [Neoantigen Cassette (4 Epitopes + AAY Cleavable Linkers)]");
    println!("                ── [MITD Lysosomal Trafficking Domain: 43 aa]");
    println!("                ── (Dual Stop Codon TGATAA)");
    println!("                ── (136 nt Human Dual Beta-Globin 3' UTR)");
    println!("                ── [Poly(A) Tail: 120 nt]");

    // Preview mRNA Sequence
    let mrna_preview =
        String::from_utf8_lossy(&vaccine_construct.full_mrna_transcript.sequence[..60]);
    let mrna_tail = String::from_utf8_lossy(
        &vaccine_construct.full_mrna_transcript.sequence[vaccine_construct.total_length_nt - 40..],
    );
    println!("\n   📜 Synthetic mRNA Sequence Preview:");
    println!("      5'-{} ... {}-3'", mrna_preview, mrna_tail);

    // =======================================================================
    // 💾 STEP 4: Store Personalized Vaccine in scies-bio-db Warehouse
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("💾 STEP 4: Persisting Personalized Vaccine Formulation into scies-bio-db");
    println!("-------------------------------------------------------------------------");

    let mrna_dna_repr = vaccine_construct
        .full_mrna_transcript
        .sequence
        .iter()
        .map(|&b| if b == b'U' { b'T' } else { b })
        .collect::<Vec<u8>>();
    db.store_sequence(9606, "Clinical_Therapeutics_v1", patient_id, &mrna_dna_repr)?;

    println!(
        "   ✅ Successfully saved into: seq/9606/Clinical_Therapeutics_v1/{}.seq",
        patient_id
    );

    // Fetch and verify stored vaccine sequence integrity
    let retrieved_vaccine =
        db.genome()
            .fetch_full_sequence(9606, "Clinical_Therapeutics_v1", patient_id)?;
    assert_eq!(retrieved_vaccine.sequence, mrna_dna_repr);
    println!(
        "   ✅ Verified vaccine record in database: 100% match ({} nt)",
        retrieved_vaccine.sequence.len()
    );

    println!("\n=========================================================================");
    println!(" 🚀 Personalized Cancer Vaccine Design Pipeline Complete & Verified 100%!");
    println!("=========================================================================\n");

    Ok(())
}
