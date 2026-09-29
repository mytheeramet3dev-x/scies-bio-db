# scies-bio-db 🗄️🧬

[![Rust](https://img.shields.io/badge/rust-2021%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-0.2.0-blue.svg)](Cargo.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![scies-series](https://img.shields.io/badge/ecosystem-scies--series-green.svg)](https://github.com/mytheeramet3dev-x)

> A high-performance hybrid biological database and sequence storage engine for the SCIES ecosystem.

`scies-bio-db` combines compact memory-mapped biological sequence storage with an embedded SQLite metadata catalog. Version 0.2 hardens the original architecture around lossless sequence preservation, corruption detection, atomic persistence, and species/assembly-scoped biological queries.

## v0.2 highlights

- `.seq v2` with `SCIESSEQ` magic, versioned 64-byte header, payload lengths and CRC32 integrity checking.
- Lossless IUPAC ambiguity preservation using a compact sidecar while retaining 2-bit A/C/G/T storage.
- Read-only compatibility with legacy `.seq v1` files.
- Atomic sequence writes using a temporary file, flush/sync and rename.
- Validation for truncated/corrupt sequence files and unsupported format versions.
- Path-component validation for assembly and chromosome identifiers.
- Assembly/species-scoped gene and variant queries.
- FASTA, GFF3 and VCF ingestion APIs with metadata integration.
- SQLite catalog for taxonomy, assemblies, reference sequences, genes, transcripts, exons, proteins and variants.

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
  + CRC32            + v1 reader   assemblies      variants/proteins
```

Sequence files remain partitioned by taxonomy and assembly:

```text
my_biodb/
├── meta.db
└── seq/
    └── 9606/
        └── GRCh38.p14/
            ├── chr1.seq
            └── chr2.seq
```

## Quick start

```rust
use scies_bio_db::BioDb;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let db = BioDb::open(Path::new("./my_biodb"))?;

    // IUPAC symbols are preserved by .seq v2.
    let dna = b"ACGTNNRYACGT";
    db.store_sequence(9606, "GRCh38", "chr21", dna)?;

    let slice = db.fetch_sequence(9606, "GRCh38", "chr21", 0, dna.len() as u64)?;
    assert_eq!(slice, dna);

    Ok(())
}
```

## Coordinate convention

SCIES internal genomic intervals use **0-based half-open coordinates**: `[start, end)`. Importers are responsible for converting source conventions such as GFF3/VCF into the internal representation.

## Storage integrity

A v2 sequence file contains a 64-byte header followed by a packed canonical sequence and a sorted ambiguity sidecar. The reader validates version, lengths, coordinate bounds and CRC32 before returning decoded data. Legacy v1 files remain readable but cannot recover ambiguous bases that were already lost by the old encoder.

See [`docs/FORMAT_SPEC.md`](docs/FORMAT_SPEC.md) for the binary layout and [`docs/API_REFERENCE.md`](docs/API_REFERENCE.md) for API details.

## Development

Run the test suite:

```bash
cargo test
```

Run Clippy with warnings denied:

```bash
cargo clippy --all-targets --all-features -- -D warnings
```

Format check:

```bash
cargo fmt --all -- --check
```

## Status

Version 0.2 focuses on correctness and integrity. Performance work such as persistent mmap caching, zero-allocation sequence views, broader connection pooling, and chromosome-scale benchmarking remains suitable for later releases after correctness invariants are locked down.

## License

MIT License. See [`LICENSE`](LICENSE).
