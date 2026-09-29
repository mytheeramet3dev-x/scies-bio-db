# Laboratory Readiness & Validation Policy

## Status

`scies-bio-db` v0.2 is engineered as a **research-grade biological data storage and retrieval component**. It is suitable for laboratory research workflows only after the installation and dataset used by that laboratory pass the validation gates below.

It is **not a medical device, diagnostic system, clinical decision-support system, or independently validated source of biological truth**. Any clinical or regulated use requires separate validation under the applicable quality and regulatory framework.

## What "lab-ready" means here

A release is considered eligible for research-laboratory deployment when all of the following hold:

1. The pinned source revision passes CI on Linux, macOS, and Windows.
2. `cargo fmt`, Clippy with warnings denied, all tests, documentation build, and release build pass.
3. Sequence round-trip tests cover canonical DNA, RNA/U, the supported IUPAC alphabet, empty sequences, odd lengths, and large random inputs.
4. Corruption tests demonstrate deterministic rejection of truncated headers, truncated payloads, invalid lengths, unsupported versions, and checksum mismatches.
5. Coordinate tests verify the documented 0-based half-open convention and boundary behavior.
6. Path traversal and malformed identifier tests pass.
7. Ingestion tests verify FASTA/GFF3/VCF parsing against fixed, versioned fixtures and expected record counts.
8. A laboratory validates its exact reference datasets against independent checksums and known-answer queries before use.
9. The exact crate version, Git commit, reference dataset versions, input hashes, and command/configuration are recorded for every analysis.
10. Backup and restore are tested on a disposable copy before the database becomes a system of record.

## Required deployment validation

For each reference assembly or annotation release, record:

- provider and release identifier;
- download date and source;
- SHA-256 (or stronger) digest of every source artifact;
- taxonomy ID and assembly accession/name;
- ingestion command and software commit;
- ingested chromosome/feature/variant counts;
- post-ingestion known-answer queries;
- operator and validation date.

Never silently replace a reference dataset in an existing analysis. Treat a changed assembly or annotation release as a new immutable dataset version.

## Data integrity model

The `.seq v2` format provides a versioned header, payload-length validation, CRC32 corruption detection, lossless IUPAC ambiguity recovery, and atomic temporary-file + sync + rename persistence. CRC32 is an integrity/error-detection checksum, **not** a cryptographic authenticity mechanism. Laboratories requiring provenance or tamper evidence should additionally maintain SHA-256 or stronger manifests outside the `.seq` payload.

SQLite metadata and sequence payloads must be backed up as one logical dataset. Do not assume a copied `meta.db` without its matching `seq/` tree is a valid backup.

## Reproducibility record

A result derived from this database should be traceable to at least:

```text
scies-bio-db version:
git commit:
OS / architecture:
Rust toolchain:
dataset provider:
dataset release:
source SHA-256:
tax_id:
assembly:
ingestion timestamp:
analysis command/config:
```

## Release rule

Do not label a revision `lab-ready` merely because it compiles. A candidate release must have green automated quality gates plus a documented validation run on representative real reference data. Failures are release blockers, not warnings.
