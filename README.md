# scies-bio-db 🗄️🧬

[![Rust](https://img.shields.io/badge/rust-2021%2B-orange.svg)](https://www.rust-lang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![scies-series](https://img.shields.io/badge/ecosystem-scies--series-green.svg)](https://github.com/mytheeramet3dev-x)

> A high-performance, hybrid biological database and sequence storage engine for all living organisms.
> Part of the **scies-series** scientific computing framework.

`scies-bio-db` bridges petabyte-scale genomic sequence storage with embedded relational metadata. Built on top of [`scies-bio-th`](https://github.com/mytheeramet3dev-x/scies-bio-th), it couples **2-bit memory-mapped sequence storage** (75% disk reduction, zero-copy random access) with an ACID-compliant **SQLite catalog** partitioned by NCBI Taxonomy ID.

---

## 🏛️ Architecture Overview

```text
                                  ┌────────────────────────┐
                                  │   BioDb Public API     │
                                  │ (Query & Ingestion Hub)│
                                  └───────────┬────────────┘
                                              │
                    ┌─────────────────────────┴─────────────────────────┐
                    ▼                                                   ▼
       ┌────────────────────────┐                          ┌────────────────────────┐
       │     SeqStore (I/O)     │                          │   SqliteStore (Meta)   │
       │  2-Bit Mmap Sequences  │                          │ Embedded SQLite Engine │
       └────────────┬───────────┘                          └────────────┬───────────┘
                    │                                                   │
     ┌──────────────┴──────────────┐                     ┌──────────────┴──────────────┐
     ▼                             ▼                     ▼                             ▼
[.seq Files]                [Zero-Copy Mmap]       [meta.db Tables]            [Relational Indices]
- 4 bases / byte            - Virtual address      - organisms                 - gene symbol index
- 8-byte LE header            translation          - chromosomes               - interval ranges
- Partitioned by TaxID      - OS page caching      - genes & transcripts       - UniProt accession
```

### Directory Hierarchy
```text
my_biodb/
├── meta.db                         # SQLite catalog with indexed schemas
└── seq/                            # Sequence payload root
    ├── 9606/                       # NCBI TaxID: Homo sapiens (Human)
    │   └── GRCh38.p14/
    │       ├── chr1.seq
    │       └── chr11_INS.seq
    ├── 2697049/                    # NCBI TaxID: SARS-CoV-2 (COVID-19)
    │   └── NC_045512.2/
    │       └── Spike_RBD.seq
    └── 2562234/                    # NCBI TaxID: Panthera spelaea (Cave Lion)
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

## 🚀 Quick Start Example

```rust
use scies_bio_db::BioDb;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Initialize or open database directory
    let db = BioDb::open(Path::new("./my_biodb"))?;

    // 2. Store DNA sequence in 2-bit compressed format
    let raw_dna = b"ATGCGATCGATCGATCGATCGATC";
    db.store_sequence(9606, "GRCh38", "chr21", raw_dna)?;

    // 3. Fetch sequence slice [0..10) via zero-copy memory mapping
    let slice = db.fetch_sequence(9606, "GRCh38", "chr21", 0, 10)?;
    println!("Decoded 10 bp: {}", String::from_utf8_lossy(&slice));

    Ok(())
}
```

---

## 🧪 Interactive Biological Simulation Demos

`scies-bio-db` includes 4 complete, end-to-end biological engineering pipelines:

### 1. 🦁 Cave Lion De-Extinction & Artificial Womb Engine
Simulates the de-extinction of the Eurasian Cave Lion (*Panthera spelaea*, TaxID: 2562234) using IVG synthetic egg design, Zona Pellucida protective matrix (ZP1-4), SCNT nuclear transfer + KDM4A histone barrier removal, and 110-day gestation in an artificial amniotic womb (Biobag) without surrogate animals:
```bash
cargo run --example cave_lion_deextinction_demo
```

### 2. 🎯 Personalized Cancer mRNA Vaccine Pipeline
Executes clinical-grade neoantigen identification: somatic mutation detection (BRAF V600E, TP53 R175H, KRAS G12D), HLA-A*02:01 MHC binding affinity scoring, polytope cassette assembly with AAY linkers, and human codon optimization:
```bash
cargo run --example cancer_vaccine_demo
```

### 3. 🥚 In Vitro Gametogenesis (IVG) & SCNT Cloning
Demonstrates stem cell differentiation into primordial germ cell-like cells (PGCLC), synthetic oocyte maturation, somatic cell nuclear transfer, and pre-implantation embryo culture:
```bash
cargo run --example ivg_scnt_demo
```

### 4. 🧬 Authentic Multi-Species DNA Ingestion & Verification
Stores and verifies genuine biological sequences from NCBI GenBank: Human Insulin (*INS*), Beta-Globin (*HBB*), *E. coli* *lacZ*, and the SARS-CoV-2 Spike Glycoprotein:
```bash
cargo run --example real_dna_demo
```

---

## 📚 Detailed Documentation

- 🏛️ **[System Architecture](docs/ARCHITECTURE.md)**: Dual-engine design, thread safety, and taxonomy partitioning.
- 💾 **[2-Bit Binary Format Spec (`.seq`)](docs/FORMAT_SPEC.md)**: Bit packing, little-endian header layout, and $O(1)$ random-access slicing equations.
- 📖 **[API Reference Guide](docs/API_REFERENCE.md)**: Comprehensive reference for `BioDb`, queries, and ingestion pipelines.
- 🦁 **[Tutorial & Pipeline Walkthroughs](docs/TUTORIAL_DEMOS.md)**: Detailed scientific rationale for the 4 simulation demos.

---

## 🧪 Testing & Verification

Run the full test suite (unit tests and doc-tests):
```bash
cargo test
```

Run clippy verification (zero warnings):
```bash
cargo clippy -- -D warnings
```

---

## 📄 License

Licensed under the [MIT License](LICENSE).
