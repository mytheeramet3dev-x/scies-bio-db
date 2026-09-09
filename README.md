# scies-bio-db 🗄️🧬

> A high-performance biological database engine for all living organisms, starting from *Homo sapiens*.
> ส่วนหนึ่งของชุดเครื่องมือวิทยาศาสตร์การคำนวณ **scies-series**

`scies-bio-db` is a hybrid biological database and sequence storage engine built on top of [`scies-bio-th`](https://github.com/mytheeramet3dev-x/scies-bio-th). It combines ultra-compact **2-bit memory-mapped sequence storage** with an embedded **SQLite relational metadata catalog**, partitioned strictly by NCBI Taxonomy ID.

---

## 🏛️ Architecture

```text
                                  scies-bio-db
                                       │
        ┌──────────────────────────────┴──────────────────────────────┐
        ▼                                                             ▼
  [Storage Engine]                                              [Catalog Layer]
  • 2-Bit Encoded DNA (.seq)                                    • SQLite Embedded (meta.db)
  • Zero-Copy Memory Mapping (memmap2)                          • Taxonomy Tree (NCBI TaxID)
  • 4 bases / byte (75% compression)                            • Chromosomes, Genes & Transcripts
  • Fast random slice queries [start..end)                      • Protein entries, Domains, Variants
```

### Directory Layout
```text
my_biodb/
├── meta.db                         # SQLite database with relational schemas
└── seq/
    ├── 9606/                       # Homo sapiens (Human)
    │   └── GRCh38/
    │       ├── chr1.seq
    │       └── chr21.seq
    ├── 2697049/                    # SARS-CoV-2 (COVID-19)
    │   └── NC_045512.2/
    │       └── Spike_RBD.seq
    └── 2562234/                    # Panthera spelaea (Cave Lion)
        └── PanSpe_Yakutia_v1.0/
            └── ExUtero_Clone_Cub_01.seq
```

---

## 📦 Package Information

- **Crate Name**: `scies-bio-db`
- **Rust Import Path**: `scies_bio_db`
- **Repository**: [https://github.com/mytheeramet3dev-x/scies-bio-db](https://github.com/mytheeramet3dev-x/scies-bio-db)
- **License**: MIT

---

## 🚀 Quick Start

```rust
use scies_bio_db::BioDb;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Open or initialize database directory
    let db = BioDb::open(Path::new("./my_biodb"))?;

    // Store a DNA sequence in 2-bit compressed format
    let raw_dna = b"ATGCGATCGATCGATCGATCGATC";
    db.store_sequence(9606, "GRCh38", "chr21", raw_dna)?;

    // Fetch sub-sequence with zero-copy random access
    let slice = db.fetch_sequence(9606, "GRCh38", "chr21", 0, 10)?;
    println!("Fetched 10 bp: {}", String::from_utf8_lossy(&slice));

    Ok(())
}
```

---

## 🧪 Interactive Production Demos

`scies-bio-db` includes 4 complete, runnable biological simulation pipelines:

### 1. 🦁 Cave Lion De-Extinction & Artificial Womb
Simulates the de-extinction of the Eurasian Cave Lion (*Panthera spelaea*, TaxID: 2562234) using IVG synthetic egg design, Zona Pellucida glycoproteins (ZP1-4), SCNT nuclear transfer + KDM4 histone barrier removal, and 110-day gestation in an artificial amniotic womb (Biobag) without surrogate mothers:
```bash
cargo run --example cave_lion_deextinction_demo
```

### 2. 🎯 Personalized Cancer mRNA Vaccine Pipeline
Executes clinical-grade neoantigen screening: somatic mutation detection (BRAF V600E, TP53 R175H, KRAS G12D), HLA-A*02:01 MHC binding affinity scoring, polytope cassette assembly with AAY linkers, and human codon optimization:
```bash
cargo run --example cancer_vaccine_demo
```

### 3. 🥚 In Vitro Gametogenesis (IVG) & SCNT
Demonstrates stem cell differentiation into primordial germ cell-like cells (PGCLC), synthetic oocyte maturation, somatic cell nuclear transfer, and pre-implantation embryo culture:
```bash
cargo run --example ivg_scnt_demo
```

### 4. 🧬 Real Biological DNA Ingestion
Stores and queries authentic sequences for Human Insulin (*INS*), Beta-Globin (*HBB*), *E. coli* *lacZ*, and the SARS-CoV-2 Spike Glycoprotein:
```bash
cargo run --example real_dna_demo
```

---

## 🧪 Testing & Quality

Run the test suite:
```bash
cargo test
```

Run clippy verification:
```bash
cargo clippy -- -D warnings
```

---

## 📄 License

Licensed under the [MIT License](LICENSE).
