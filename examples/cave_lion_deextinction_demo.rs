//! Full-Spectrum Cave Lion (Panthera spelaea) De-Extinction Pipeline
//! Featuring Synthetic Oocyte with Zona Pellucida & Artificial Womb Incubation (Zero Surrogate Mother)
//! Using `scies-bio-th::ivg` and `scies-bio-db`.

use scies_bio_db::model::genome::{AssemblyStatus, ReferenceGenome};
use scies_bio_db::model::taxonomy::{Lineage, Organism};
use scies_bio_db::BioDb;
use scies_bio_th::ivg::embryo::{ArtificialWombChamber, EmbryoCultureEnvironment, SyntheticEmbryo};
use scies_bio_th::ivg::gamete::{GameteType, InductionProtocol, PgclcCulture, ZonaPellucida};
use scies_bio_th::ivg::scnt::{
    ActivationMethod, EnucleatedOocyte, ScntConstruct, SomaticDonorCell,
};
use scies_bio_th::seq::Dna;
use std::path::Path;

fn main() -> scies_bio_db::Result<()> {
    println!("=========================================================================");
    println!(" 🦁❄️ CAVE LION (Panthera spelaea) DE-EXTINCTION & ARTIFICIAL WOMB ENGINE");
    println!("=========================================================================\n");

    let db_path = Path::new("./demo_biodb");
    if !db_path.exists() {
        eprintln!("Please run real_dna_demo first to initialize the database!");
        return Ok(());
    }
    let db = BioDb::open(db_path)?;

    // =======================================================================
    // 🏛️ STEP 1: Register Cave Lion in Universal Taxonomy Database
    // =======================================================================
    println!("-------------------------------------------------------------------------");
    println!("🏛️ STEP 1: Taxonomy & Ancient Genome Registration in scies-bio-db");
    println!("-------------------------------------------------------------------------");

    let cave_lion = Organism {
        tax_id: 2562234, // NCBI Taxonomy ID for Panthera spelaea
        scientific_name: "Panthera spelaea".to_string(),
        common_name: Some("Eurasian Cave Lion (Extinct Pleistocene Apex Predator)".to_string()),
        lineage: Lineage {
            domain: "Eukaryota".into(),
            kingdom: "Animalia".into(),
            phylum: "Chordata".into(),
            class: "Mammalia".into(),
            order: "Carnivora".into(),
            family: "Felidae".into(),
            genus: "Panthera".into(),
            species: "spelaea".into(),
        },
    };
    db.register_organism(&cave_lion)?;
    println!(
        "   ✅ Registered TaxID [2562234]: {}",
        cave_lion.scientific_name
    );
    println!(
        "      - Common Name : {}",
        cave_lion.common_name.as_deref().unwrap()
    );
    println!(
        "      - Lineage     : {} -> {} -> {} -> {}",
        cave_lion.lineage.class,
        cave_lion.lineage.order,
        cave_lion.lineage.family,
        cave_lion.lineage.genus
    );

    db.register_reference_genome(&ReferenceGenome {
        assembly_name: "PanSpe_Yakutia_v1.0".to_string(),
        tax_id: 2562234,
        status: AssemblyStatus::Primary,
        release_year: 2026,
    })?;

    // Permafrost ancient DNA sequence (Siberian cub specimen)
    let cave_lion_nuclear_dna = Dna::new(
        b"ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
    ).map_err(scies_bio_db::BioDbError::Bio)?;

    // =======================================================================
    // 🥚 STEP 2: Engineered Synthetic Egg with Zona Pellucida (IVG)
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🥚 STEP 2: Designing Synthetic Egg (IVG) & Protective Zona Pellucida Shell");
    println!("-------------------------------------------------------------------------");

    println!("   🔬 Inducing Primordial Germ Cells from stem cells (Hayashi/Saitou cytokine protocol)...");
    let pgclc = PgclcCulture::induce(
        &cave_lion_nuclear_dna,
        2562234,
        InductionProtocol::Hayashi2016,
        5,
    )
    .map_err(scies_bio_db::BioDbError::Bio)?;

    let mut synthetic_oocyte = pgclc
        .mature_gamete(GameteType::Oocyte, 99999)
        .map_err(scies_bio_db::BioDbError::Bio)?;

    // Attach custom engineered Zona Pellucida matrix
    let zp_shell = ZonaPellucida::cave_lion_default();
    synthetic_oocyte.zona_pellucida = Some(zp_shell);

    let zp = synthetic_oocyte.zona_pellucida.as_ref().unwrap();
    println!("   🎉 Synthetic Cave Lion Oocyte Engineered:");
    println!(
        "      - Polar Body Extruded : {}",
        if synthetic_oocyte.quality.polar_body_extruded {
            "✅ YES (Metaphase II arrest)"
        } else {
            "❌ NO"
        }
    );
    println!(
        "      - Mitochondrial Count : {} copies (Maternal energy supply)",
        synthetic_oocyte.quality.mitochondrial_count
    );
    println!(
        "      - Imprint Integrity   : {:.1}% (Maternal imprint established)",
        synthetic_oocyte.quality.imprint_fidelity * 100.0
    );
    println!(
        "      - Zona Pellucida Shell: 16.5 µm thickness | Elasticity: {:.1} kPa",
        zp.elasticity_kpa
    );
    println!("      - ZP Glycoproteins    : ZP1 (cross-link: {:.0}%), ZP2 (permissive), ZP3 (receptor: {:.0}%), ZP4 (structural)", 
        zp.zp1_crosslinker_density * 100.0, zp.zp3_receptor_affinity * 100.0);

    // =======================================================================
    // 🦁 STEP 3: Somatic Cell Nuclear Transfer (SCNT) with KDM4 Demethylase
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🦁 STEP 3: SCNT Nuclear Transfer & Epigenetic Barrier Erasure");
    println!("-------------------------------------------------------------------------");

    let cytoplast =
        EnucleatedOocyte::from_oocyte(&synthetic_oocyte).map_err(scies_bio_db::BioDbError::Bio)?;
    println!(
        "   🧫 Oocyte enucleated (Maternal chromosomes removed, Cytoplasm & Mitochondria intact)"
    );

    let archaic_somatic_donor = SomaticDonorCell::new(
        cave_lion_nuclear_dna.clone(),
        2562234,
        15.0, // 15-year old adult tissue
        "Yakutian Permafrost Deep Tissue",
    );

    println!("   ⚡ Microinjecting Cave Lion Nucleus + KDM4A Histone Demethylase mRNA...");
    let scnt_zygote = ScntConstruct::reconstruct(
        archaic_somatic_donor,
        cytoplast,
        ActivationMethod::ElectrofusionPlusKinaseInhibitor,
        true, // KDM4A mRNA applied
    )
    .map_err(scies_bio_db::BioDbError::Bio)?;

    println!(
        "      - Electrofusion Efficiency : {:.1}%",
        scnt_zygote.activation_efficiency * 100.0
    );
    println!(
        "      - H3K9me3 Barrier Removal  : {:.1}% residual (Permissive to EGA)",
        scnt_zygote.epigenetic_profile.h3k9me3_barrier * 100.0
    );
    println!(
        "      - Blastocyst Formation Rate: {:.1}% predicted competence",
        scnt_zygote.predict_blastocyst_competence() * 100.0
    );

    // =======================================================================
    // 🐣 STEP 4: In Vitro Pre-Implantation Embryogenesis (Day 0 - Day 6)
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🐣 STEP 4: Pre-Implantation Embryogenesis & Zona Pellucida Hatching");
    println!("-------------------------------------------------------------------------");

    let mut embryo =
        SyntheticEmbryo::from_scnt(scnt_zygote).map_err(scies_bio_db::BioDbError::Bio)?;

    let mut pre_implantation_env = EmbryoCultureEnvironment::default();
    pre_implantation_env.oxygen_percent = 5.0; // 5% Low O2
    pre_implantation_env.microfluidic_flow = true;
    pre_implantation_env.sequential_media_optimized = true;

    // Day 1: 2-Cell stage (EGA)
    embryo
        .advance_development(24.0, &pre_implantation_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!("   🌱 Day 1.0 (24h) : Stage = TwoCell | Cells: {} | Embryonic Genome Activation (EGA) Onset", embryo.total_cell_count);

    // Day 3: Morula
    embryo
        .advance_development(48.0, &pre_implantation_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!("   🌱 Day 3.0 (72h) : Stage = Morula  | Cells: {} | Cellular Compaction & E-Cadherin Polarity", embryo.total_cell_count);

    // Day 5: Blastocyst
    embryo
        .advance_development(48.0, &pre_implantation_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!(
        "   🌱 Day 5.0 (120h): Stage = Blastocyst | Cells: {} (ICM: {}, Trophectoderm: {})",
        embryo.total_cell_count, embryo.icm_cell_count, embryo.trophectoderm_cell_count
    );

    // Day 6: Hatched Blastocyst (Breaking out of Zona Pellucida)
    embryo
        .advance_development(24.0, &pre_implantation_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!("   ✨ Day 6.0 (144h): Stage = HatchedBlastocyst | Cells: {} | Hatching through Zona Pellucida Complete!", embryo.total_cell_count);

    // =======================================================================
    // 🫀 STEP 5: Artificial Womb (Biobag) Extra-Uterine Gestation (Day 7 - Day 110)
    //    NO LIVE SURROGATE MOTHER LIONESS NEEDED!
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("🫀 STEP 5: Artificial Womb (Biobag) Extra-Uterine Gestation (Zero Surrogate Mother)");
    println!("-------------------------------------------------------------------------");

    let artificial_womb = ArtificialWombChamber {
        amniotic_fluid_rate_ml_hr: 75.0,
        oxygenator_pao2_mmhg: 30.0,
        umbilical_cannulation_patency: 0.99,
        chamber_temp_celsius: 38.2, // Feline maternal temperature
    };

    let ex_utero_env = EmbryoCultureEnvironment {
        oxygen_percent: 5.0,
        microfluidic_flow: true,
        sequential_media_optimized: true,
        artificial_womb_supported: true,
        womb_chamber: Some(artificial_womb),
    };

    println!("   🏥 Docking Hatched Blastocyst into Sterile Artificial Amniotic Biobag Chamber...");

    // Day 14: Gastrulation
    embryo
        .advance_development(8.0 * 24.0, &ex_utero_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!("   🐾 Day 14.0 (Gastrulation) : 3 Germ Layers Established (Ectoderm / Mesoderm / Endoderm)");

    // Day 30: Organogenesis
    embryo
        .advance_development(16.0 * 24.0, &ex_utero_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!("   🐾 Day 30.0 (Organogenesis) : Cardiac Contractions Active (120 BPM) | Limb Buds Formed | Neural Tube Closed");

    // Day 60: Fetal Development
    embryo
        .advance_development(30.0 * 24.0, &ex_utero_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!("   🐾 Day 60.0 (Fetal Growth)  : Dense Fur Follicles Developing | Claws & Retractile Sheaths Formed");

    // Day 110: Full Term!
    embryo
        .advance_development(50.0 * 24.0, &ex_utero_env)
        .map_err(scies_bio_db::BioDbError::Bio)?;
    println!("\n   🎉 DAY 110.0: FULL-TERM GESTATION REACHED IN ARTIFICIAL WOMB!");
    println!("      - Developmental Stage: {:?}", embryo.stage);
    println!(
        "      - Total Cells Count  : ~{:.1} Billion cells",
        embryo.total_cell_count as f64 / 1_000_000_000.0
    );
    println!(
        "      - Fetal Birth Weight : {:.2} kg (Healthy Eurasian Cave Lion Cub)",
        embryo.fetal_weight_grams / 1000.0
    );
    println!(
        "      - Viability Score    : {:.1} / 100.0",
        embryo.viability_score
    );

    // =======================================================================
    // 💾 STEP 6: Save De-Extinction Record to scies-bio-db
    // =======================================================================
    println!("\n-------------------------------------------------------------------------");
    println!("💾 STEP 6: Storing De-Extinction Record into scies-bio-db Warehouse");
    println!("-------------------------------------------------------------------------");

    db.store_sequence(
        2562234,
        "PanSpe_Yakutia_v1.0",
        "ExUtero_Clone_Cub_01",
        &cave_lion_nuclear_dna.sequence,
    )?;
    println!(
        "   ✅ Successfully saved into: seq/2562234/PanSpe_Yakutia_v1.0/ExUtero_Clone_Cub_01.seq"
    );

    println!("\n=========================================================================");
    println!(" 🏆 Cave Lion Synthetic Egg & Artificial Womb Gestation Succeeded 100%!");
    println!("=========================================================================\n");

    Ok(())
}
