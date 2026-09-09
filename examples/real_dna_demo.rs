//! Demonstration of storing and querying real biological DNA sequences
//! across multiple organisms in `scies-bio-db`.

use scies_bio_db::model::genome::{AssemblyStatus, ReferenceGenome};
use scies_bio_db::model::taxonomy::{Lineage, Organism};
use scies_bio_db::BioDb;
use std::fs;
use std::path::Path;

fn main() -> scies_bio_db::Result<()> {
    println!("===============================================================");
    println!(" 🧬 scies-bio-db: Real Biological DNA Storage Demonstration");
    println!("===============================================================\n");

    let db_path = Path::new("./demo_biodb");
    // Clean up previous demo run if any
    let _ = fs::remove_dir_all(db_path);

    // 1. Initialize the database
    println!("1️⃣  Opening biological database at: {:?}", db_path);
    let db = BioDb::open(db_path)?;

    // -----------------------------------------------------------------------
    // 2. Register Taxonomy for Multiple Kingdoms & Domains
    // -----------------------------------------------------------------------
    println!("\n2️⃣  Registering Organisms in Taxonomy...");

    // Human (Homo sapiens)
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
    db.register_organism(&human)?;
    println!("   ✅ [TaxID: 9606] {}", human.scientific_name);

    // E. coli (Escherichia coli K-12)
    let ecoli = Organism {
        tax_id: 562,
        scientific_name: "Escherichia coli".to_string(),
        common_name: Some("E. coli (Bacteria)".to_string()),
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
    db.register_organism(&ecoli)?;
    println!("   ✅ [TaxID: 562]  {}", ecoli.scientific_name);

    // SARS-CoV-2 (Virus)
    let virus = Organism {
        tax_id: 2697049,
        scientific_name: "Severe acute respiratory syndrome coronavirus 2".to_string(),
        common_name: Some("SARS-CoV-2 (Virus)".to_string()),
        lineage: Lineage {
            domain: "Viruses".into(),
            kingdom: "Orthornavirae".into(),
            phylum: "Pisuviricota".into(),
            class: "Pisoniviricetes".into(),
            order: "Nidovirales".into(),
            family: "Coronaviridae".into(),
            genus: "Betacoronavirus".into(),
            species: "SARS-CoV-2".into(),
        },
    };
    db.register_organism(&virus)?;
    println!("   ✅ [TaxID: 2697049] {}", virus.scientific_name);

    // -----------------------------------------------------------------------
    // 3. Register Assemblies
    // -----------------------------------------------------------------------
    println!("\n3️⃣  Registering Reference Assemblies...");
    db.register_reference_genome(&ReferenceGenome {
        assembly_name: "GRCh38.p14".to_string(),
        tax_id: 9606,
        status: AssemblyStatus::Primary,
        release_year: 2013,
    })?;

    db.register_reference_genome(&ReferenceGenome {
        assembly_name: "ASM584v2".to_string(),
        tax_id: 562,
        status: AssemblyStatus::Primary,
        release_year: 2014,
    })?;

    db.register_reference_genome(&ReferenceGenome {
        assembly_name: "NC_045512.2".to_string(),
        tax_id: 2697049,
        status: AssemblyStatus::Primary,
        release_year: 2020,
    })?;

    // -----------------------------------------------------------------------
    // 4. Create FASTA files with Real DNA Sequences from NCBI
    // -----------------------------------------------------------------------
    println!("\n4️⃣  Preparing Real DNA Sequences (NCBI/RefSeq)...");

    let tmp_dir = Path::new("./demo_biodb/raw_fasta");
    fs::create_dir_all(tmp_dir)?;

    // A. Human Insulin (INS) & Hemoglobin Beta (HBB) coding sequences
    let human_fasta = tmp_dir.join("human_genes.fa");
    fs::write(
        &human_fasta,
        b">chr11_INS Homo sapiens insulin (INS) CDS (NM_000207.3)\n\
        ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGAC\n\
        CCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTAC\n\
        CTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGAC\n\
        CTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTG\n\
        GCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGC\n\
        TCCCTCTACCAGCTGGAGAACTACTGCAACTAG\n\
        >chr11_HBB Homo sapiens hemoglobin subunit beta (HBB) CDS (NM_000518.5)\n\
        ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC\n\
        GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAG\n\
        AGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAG\n\
        GTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGAC\n\
        AACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT\n\
        CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC\n\
        AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAAT\n\
        GCCCTGGCCCACAAGTATCACTAA\n",
    )?;

    // B. E. coli lacZ 5' region (NC_000913.3)
    let ecoli_fasta = tmp_dir.join("ecoli_lacZ.fa");
    fs::write(
        &ecoli_fasta,
        b">lacZ Escherichia coli str. K-12 substr. MG1655 beta-galactosidase\n\
        ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCT\n\
        GGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGC\n\
        GAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC\n\
        TTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCT\n\
        GAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATC\n\
        TACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCG\n\
        ACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACG\n\
        CGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTC\n\
        GGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCC\n\
        GGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGGCAGTTATCTGGAAGAT\n\
        CAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACT\n\
        ACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTA\n\
        CTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCT\n\
        TTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATC\n\
        GATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAACGTCGAAAACCCGAAA\n\
        CTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGAC\n\
        GGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGAGGTGCGGATTGAAAAT\n\
        GGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGTTAACCGTCACGAGCAT\n\
        CATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATG\n\
        AAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAACCATCCGCTGTGGTAC\n\
        ACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGC\n\
        ATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACCGGCGATGAGCGAACGC\n\
        GTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGATCATCTGGTCGCTGGGG\n\
        AATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGAT\n\
        CCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTATT\n\
        TGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCC\n\
        ATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGATCCTTTGCGAATACGCC\n\
        CACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCAGGCGTTTCGTCAGTAT\n\
        CCCCGTTTACAGGGCGGCTTCGTCTGGGACTGGGTGGATCAGTCGCTGATTAAATATGAT\n\
        GAAAACGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGATACGCCGAACGATCGC\n\
        CAGTTCTGTATGAACGGTCTGGTCTTTGCCGACCGCACGCCGCATCCAGCGCTGACGGAA\n\
        GCAAAACACCAGCAGCAGTTTTTCCAGTTCCGTTTATCCGGGCAAACCATCGAAGTGACC\n\
        AGCGAATACCTGTTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGAT\n\
        GGTAAGCCGCTGGCAAGCGGTGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAGTTG\n\
        ATTGAACTGCCTGAACTACCGCAGCCGGAGAGCGCCGGGCAACTCTGGCTCACAGTACGC\n\
        GTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACATCAGCGCCTGGCAGCAG\n\
        TGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCGCGTCCCACGCCATCCCGCAT\n\
        CTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAAGCGTTGGCAATTTAAC\n\
        CGCCAGTCAGGCTTTCTTTCACAGATGTGGATTGGCGATAAAAAACAACTGCTGACGCCG\n\
        CTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACATTGGCGTAAGTGAAGCGACC\n\
        CGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGGCCATTACCAGGCCGAA\n\
        GCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGATGCGGTGCTGATTACGACCGCT\n\
        CACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAAAACCTACCGGATTGAT\n\
        GGTAGTGGTCAAATGGCGATTACCGTTGATGTTGAAGTGGCGAGCGATACACCGCATCCG\n\
        GCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCGGGTAAACTGGCTCGGA\n\
        TTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTGTTTTGACCGCTGGGAT\n\
        CTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGAGCGAAAACGGTCTGCGCTGC\n\
        GGGACGCGCGAATTGAATTATGGCCCACACCAGTGGCGCGGCGACTTCCAGTTCAACATC\n\
        AGCCGCTACAGTCAACAGCAACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAA\n\
        GAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGGTGGCGACGACTCCTGG\n\
        AGCCCGTCAGTATCGGCGGAATTCCAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTC\n\
        TGGTGTCAAAAATAA\n",
    )?;

    // C. SARS-CoV-2 Spike Receptor Binding Domain (RBD)
    let sars_fasta = tmp_dir.join("sars_cov2_spike_rbd.fa");
    fs::write(
        &sars_fasta,
        b">Spike_RBD SARS-CoV-2 Spike Glycoprotein RBD Region (NC_045512.2)\n\
        AACATCACAAATCTTTGTCCTTTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTT\n\
        TATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAAT\n\
        TCCGCATCATTTTCCACTTTTAAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTC\n\
        TGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATC\n\
        GCTCCAGGGCAAACTGGAAAGATTGCTGATTATAATTATAAATTACCAGATGATTTTACA\n\
        GGCTGCGTTATAGCTTGGAATTCTAACAATCTTGATTCTAAGGTTGGTGGTAATTATAAT\n\
        TACCTGTATAGATTGTTTAGGAAGTCTAATCTCAAACCTTTTGAGAGAGATATTTCAACT\n\
        GAAATCTATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTT\n\
        CCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGTA\n\
        GTAGTACTTTCTTTTGAACTTCTACATGCACCAGCAACTGTTTGTGGACCTAAAAAGTCT\n\
        ACTAATTTGGTTAAAAACAAATGTGTCAATTTCAACTTCAATGGTTTAACAGGCACAGGT\n\
        GTTCTTACTGAGTCTAACAAAAAGTTTCTGCCTTTCCAACAATTTGGCAGAGACATTGCT\n\
        GACACTACTGATGCTGTCCGTGATCCACAGACACTTGAGATTCTTGACATTACACCATGT\n\
        TCTTTTGGTGGTGTCAGTGTTATAACACCAGGAACAAATACTTCTAACCAGGTTGCTGTT\n\
        CTTTATCAGGATGTTAACTGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACT\n\
        CCTACTTGGCGTGTTTATTCTACAGGTTCTAATGTTTTTCAAACACGTGCAGGCTGTTTA\n\
        ATAGGGGCTGAACATGTCAACAACTCATATGAGTGTGACATACCCATTGGTGCAGGTATA\n\
        TGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAA\n\
        TCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTACTCTAATAAC\n\
        TCTATTGCCATACCCACAAATTTTACTATTAGTGTTACCACAGAAATTCTACCAGTGTCT\n\
        ATGACCAAGACATCAGTAGATTGTACAATGTACATTTGTGGTGATTCAACTGAATGCAGC\n\
        AATCTTTTGTTGCAATATGGCAGTTTTTGTACACAATTAAACCGTGCTTTAACTGGAATA\n\
        GCTGTTGAACAAGACAAAAACACCCAAGAAGTTTTTGCACAAGTCAAACAAATTTACAAA\n\
        ACACCACCAATTAAAGATTTTGGTGGTTTTAATTTTTCACAAATATTACCAGATCCATCA\n\
        AAACCAAGCAAGAGGTCATTTATTGAAGATCTACTTTTCAACAAAGTGACACTTGCAGAT\n\
        GCTGGCTTCATCAAACAATATGGTGATTGCCTTGGTGATATTGCTGCTAGAGACCTCATTT\n\
        GTGCACAAAAGTTTAACGGCCTTACTGTTTTGCCACCTTTGCTCACAGATGAAATGATTG\n\
        CTCAATACACTTCTGCACTGTTAGCGGGTACAATCACTTCTGGTTGGACCTTTGGTGCAG\n\
        GTGCTGCATTACAAATACCATTTGCTATGCAAATGGCTTATAGGTTTAATGGTATTGGAG\n\
        TTACACAGAATGTTCTCTATGAGAACCAAAAATTGATTGCCAACCAATTTAATAGTGCTA\n\
        TTGGCAAAATTCAAGACTCACTTTCTTCCACAGCAAGTGCACTTGGAAAACTTCAAGATG\n\
        TGGTCAACCAAAATGCACAAGCTTTAAACACGCTTGTTAAACAACTTAGCTCCAATTTTG\n\
        GTGCAATTTCAAGTGTTTTAAATGATATCCTTTCACGTCTTGACAAAGTTGAGGCTGAAG\n\
        TGCAAATTGATAGGTTGATCACAGGCAGACTTCAAAGTTTGCAGACATATGTGACTCAAC\n\
        AATTAATTAGAGCTGCAGAAATCAGAGCTTCTGCTAATCTTGCTGCTACTAAAATGTCAG\n\
        AGTGTGTACTTGGACAATCAAAAAGAGTTGATTTTTGTGGAAAGGGCTATCATCTTATGT\n\
        CCTTCCCTCAGTCAGCACCTCATGGTGTAGTCTTCTTGCATGTGACTTATGTCCCTGCAC\n\
        AAGAAAAGAACTTCACAACTGCTCCTGCCATTTGTCATGATGGAAAAGCACACTTTCCTC\n\
        GTGAAGGTGTCTTTGTTTCAAATGGCACACACTGGTTTGTAACACAAAGGAATTTTTATG\n\
        AACCACAAATCATTACTACAGACAACACATTTGTGTCTGGTAACTGTGATGTTGTAATAG\n\
        GAATTGTCAACAACACAGTTTATGATCCTTTGCAACCTGAATTAGACTCATTCAAGGAGG\n\
        AGTTAGATAAATATTTTAAGAATCATACATCACCAGATGTTGATTTAGGTGACATCTCTG\n\
        GCATTAATGCTTCAGTTGTAAACATTCAAAAAGAAATTGACCGCCTCAATGAGGTCGCCA\n\
        AAAATTTAAATGAATCACTCATTGACCTTCAAGAATTGGGAAAATATGAGCAATATATTA\n\
        AATGGCCTTGGTACATTTGGCTAGGTTTTATAGCTGGCTTGATTGCCATAGTAATGGTGA\n\
        CAATTATGCTTTGCTGTATGACCAGTTGCTGTAGTTGTCTCAAGGGCTGTTGTTCTTGTC\n\
        CATCAGGCTGCTGTGATTTTGATGAAGATGACTCTGAGCCAGTTCTCAAGGGTGTCAAAT\n\
        TACATTACACATAA\n",
    )?;

    // -----------------------------------------------------------------------
    // 5. Ingest Real Sequences into scies-bio-db
    // -----------------------------------------------------------------------
    println!("\n5️⃣  Ingesting FASTA records into Database...");

    let human_ingested = db.ingest_fasta(&human_fasta, 9606, "GRCh38.p14")?;
    println!(
        "   📥 Human: ingested {} chromosome/genes: {:?}",
        human_ingested.len(),
        human_ingested
    );

    let ecoli_ingested = db.ingest_fasta(&ecoli_fasta, 562, "ASM584v2")?;
    println!(
        "   📥 E. coli: ingested {} chromosome/genes: {:?}",
        ecoli_ingested.len(),
        ecoli_ingested
    );

    let sars_ingested = db.ingest_fasta(&sars_fasta, 2697049, "NC_045512.2")?;
    println!(
        "   📥 SARS-CoV-2: ingested {} chromosome/genes: {:?}",
        sars_ingested.len(),
        sars_ingested
    );

    // -----------------------------------------------------------------------
    // 6. Query and Inspect Sequences Back from the Database
    // -----------------------------------------------------------------------
    println!("\n6️⃣  Querying Real DNA from Storage & Computing Properties...");

    // Query Human Insulin
    let ins_dna = db
        .genome()
        .fetch_full_sequence(9606, "GRCh38.p14", "chr11_INS")?;
    println!("\n   📌 [Human Insulin CDS (INS)]");
    println!("      - Length: {} bp", ins_dna.sequence.len());
    println!("      - GC Content: {:.2}%", ins_dna.gc_content());
    println!(
        "      - First 30 bp: {}",
        String::from_utf8_lossy(&ins_dna.sequence[..30])
    );
    println!(
        "      - Reverse Complement (first 30 bp): {}",
        String::from_utf8_lossy(&ins_dna.reverse_complement().sequence[..30])
    );

    // Query Human Hemoglobin Beta
    let hbb_dna = db
        .genome()
        .fetch_full_sequence(9606, "GRCh38.p14", "chr11_HBB")?;
    println!("\n   📌 [Human Hemoglobin Beta (HBB)]");
    println!("      - Length: {} bp", hbb_dna.sequence.len());
    println!("      - GC Content: {:.2}%", hbb_dna.gc_content());
    println!(
        "      - First 30 bp: {}",
        String::from_utf8_lossy(&hbb_dna.sequence[..30])
    );

    // Query E. coli lacZ
    let lacz_dna = db.genome().fetch_full_sequence(562, "ASM584v2", "lacZ")?;
    println!("\n   📌 [E. coli Beta-Galactosidase (lacZ)]");
    println!("      - Length: {} bp", lacz_dna.sequence.len());
    println!("      - GC Content: {:.2}%", lacz_dna.gc_content());
    println!(
        "      - First 30 bp: {}",
        String::from_utf8_lossy(&lacz_dna.sequence[..30])
    );

    // Query SARS-CoV-2 Spike RBD
    let spike_dna = db
        .genome()
        .fetch_full_sequence(2697049, "NC_045512.2", "Spike_RBD")?;
    println!("\n   📌 [SARS-CoV-2 Spike Glycoprotein gene]");
    println!("      - Length: {} bp", spike_dna.sequence.len());
    println!("      - GC Content: {:.2}%", spike_dna.gc_content());
    println!(
        "      - First 30 bp: {}",
        String::from_utf8_lossy(&spike_dna.sequence[..30])
    );

    // -----------------------------------------------------------------------
    // 7. Check 2-bit Compression Efficiency on Disk
    // -----------------------------------------------------------------------
    println!("\n7️⃣  2-bit Storage Compression Efficiency:");
    let total_bases = ins_dna.sequence.len()
        + hbb_dna.sequence.len()
        + lacz_dna.sequence.len()
        + spike_dna.sequence.len();

    let ins_file = db_path.join("seq/9606/GRCh38.p14/chr11_INS.seq");
    let lacz_file = db_path.join("seq/562/ASM584v2/lacZ.seq");
    let spike_file = db_path.join("seq/2697049/NC_045512.2/Spike_RBD.seq");

    let compressed_bytes = fs::metadata(&ins_file)?.len()
        + fs::metadata(db_path.join("seq/9606/GRCh38.p14/chr11_HBB.seq"))?.len()
        + fs::metadata(&lacz_file)?.len()
        + fs::metadata(&spike_file)?.len();

    println!(
        "      - Total raw bases stored: {} bp ({} bytes raw ASCII text)",
        total_bases, total_bases
    );
    println!(
        "      - Total compressed size on disk: {} bytes (including 8-byte headers)",
        compressed_bytes
    );
    println!(
        "      - Compression Ratio: ~{:.1}x smaller than raw text file!",
        total_bases as f64 / compressed_bytes as f64
    );

    println!("\n===============================================================");
    println!(" 🎉 Real DNA Storage & Multi-Species Query Succeeded 100%!");
    println!("===============================================================\n");

    Ok(())
}
