# scies-bio-db 🗄️🧬

[![Rust](https://img.shields.io/badge/rust-2021%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-0.2.0-blue.svg)](Cargo.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![scies-series](https://img.shields.io/badge/ecosystem-scies--series-green.svg)](https://github.com/mytheeramet3dev-x)

> A high-performance hybrid biological database and sequence storage engine for the SCIES ecosystem.

`scies-bio-db` combines compact memory-mapped 2-bit biological sequence storage with an embedded SQLite metadata catalog. Version 0.2 hardens the architecture around lossless sequence preservation, SHA-256 cryptographic provenance, annotation ingestion, and curated longevity intervention loci for downstream simulation in `scies-bio-th` and `Humanz`.

## Highlights

- **`.seq v2` Storage Format**: `SCIESSEQ` magic, versioned 64-byte header, CRC32 integrity checking, and lossless IUPAC sidecars preserving degenerate symbols and Uracil.
- **Cryptographic Provenance**: Every ingested reference sequence and assembly computes and records its SHA-256 hash in SQLite (`BioDb::verify_sequence_integrity`).
- **Annotation Pipeline**: Ingests GFF3 and GTF annotations with diagnostics (`Gff3IngestReport`), attribute resolution (`gene_name`, `Name`, `gene`, `symbol`), and malformed coordinate validation.
- **Longevity Target Catalog**: Curated database table seeded with 10 essential human longevity genes: `TERT`, `TERC`, `SIRT6`, `FOXO3`, `TP53`, `LMNA`, `CDKN2A`, `ATM`, `PARP1`, `KLOTHO`.
- **Typed LongevityLocus Extraction**: Retrieves full locus structures (`LongevityLocus`) combining `GeneFeature`, `scies_bio_th::Dna` sequence slices, window coordinates, and provenance hashes.
- **Downstream Ready**: Feeds directly into `scies-bio-th` CRISPR guide design and `Humanz` genome architecture simulation.
- **Path Isolation & Atomic Writes**: Atomic temporary file renames and strict traversal protection across assembly and sequence paths.

## Architecture

```text
                         BioDb public API
                                │
                 ┌─────────────┴─────────────┐
                 ▼                           ▼
          SeqStore (.seq v2)          SqliteStore
          mmap-backed reads           metadata/catalog
                 │                           │
       ┌─────────┴─────────┐       ┌─────────┴──────────┐
       ▼                   ▼       ▼                    ▼
  2-bit payload      IUPAC sidecar taxonomy       annotations
  + CRC32            + v1 reader   assemblies      longevity_targets
                                   provenance      genes / transcripts
```

Sequence files are partitioned by taxonomy and assembly:

```text
my_biodb/
├── meta.db
└── seq/
    └── 9606/
        └── T2T-CHM13v2.0/
            ├── chr5.seq
            ├── chr17.seq
            └── chr19.seq
```

## Quick Start

### 1. Ingest Reference Sequences and Annotations

```rust
use scies_bio_db::BioDb;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let db = BioDb::open(Path::new("./my_biodb"))?;

    // Ingest FASTA with cryptographic SHA-256 provenance
    db.ingest_fasta_with_uri(
        Path::new("chr5.fasta"),
        9606,
        "T2T-CHM13v2.0",
        Some("https://api.ncbi.nlm.nih.gov/genomes/T2T-CHM13v2.0/chr5.fa"),
    )?;

    // Ingest GFF3 annotations
    let report = db.ingest_gff3_report(Path::new("genes.gff3"), 9606, "T2T-CHM13v2.0")?;
    println!("Ingested {} genes", report.genes_count);

    // Query gene by symbol
    let tert_gene = db.get_gene_by_symbol(9606, "T2T-CHM13v2.0", "TERT")?;
    println!("TERT location: {}:{}-{}", tert_gene.chrom, tert_gene.start_bp, tert_gene.end_bp);

    // Query longevity targets
    let targets = db.list_longevity_targets(9606)?;
    for t in targets {
        println!("{}: prio {} [{}]", t.symbol, t.priority, t.category);
    }

    // Extract typed LongevityLocus
    let locus = db.get_longevity_locus(9606, "T2T-CHM13v2.0", "TERT")?;
    println!("Extracted {} bp with provenance hash {}", locus.len(), locus.provenance_sha256);

    Ok(())
}
```

### 2. Run the Longevity Query Example

```bash
cargo run --release --example ingest_and_query_longevity
```

## Coordinate Convention

SCIES internal genomic intervals use **0-based half-open coordinates**: `[start, end)`. Importers automatically convert source 1-based inclusive conventions from GFF3 and VCF into the internal 0-based representation.

## Storage Integrity & Verification

A `.seq v2` file contains a 64-byte header followed by packed canonical 2-bit nucleotides and a sorted ambiguity sidecar. The reader verifies version, lengths, coordinate bounds and CRC32 before decoding. Additionally, `db.verify_sequence_integrity(tax_id, assembly, seq_id)` validates the on-disk sequence payload against its SQLite SHA-256 catalog hash.

## Development

Run the test suite:

```bash
cargo test --release
```

Check all targets:

```bash
cargo check --release --all-targets
```

## License

MIT License. See [`LICENSE`](LICENSE).
