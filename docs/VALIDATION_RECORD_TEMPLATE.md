# scies-bio-db Validation Record

## Software identity

- scies-bio-db version:
- Git commit SHA:
- Rust version (`rustc --version --verbose`):
- OS / architecture:
- Validation date:
- Validator:

## Source dataset identity

- Provider:
- Dataset / reference name:
- Release / accession:
- NCBI Taxonomy ID:
- Assembly:
- Source retrieval date:
- Source artifact SHA-256:

## Automated gates

- [ ] `cargo fmt --all -- --check`
- [ ] `cargo clippy --all-targets --all-features -- -D warnings`
- [ ] `cargo test --all-features --locked`
- [ ] `RUSTDOCFLAGS="-Dwarnings" cargo doc --no-deps --all-features`
- [ ] `cargo build --release --all-features --locked`

CI run URL / ID:

## Ingestion validation

- Input format(s): FASTA / GFF3 / VCF / other
- Exact ingestion command/configuration:
- Expected chromosome count:
- Observed chromosome count:
- Expected feature count:
- Observed feature count:
- Expected variant count:
- Observed variant count:
- Errors/warnings:

## Known-answer checks

Record queries whose expected answers were independently established.

| Test | Expected | Observed | Pass |
|---|---|---|---|
| sequence boundary query | | | |
| ambiguity/IUPAC query | | | |
| gene lookup | | | |
| regional lookup | | | |
| variant lookup | | | |

## Integrity / recovery

- [ ] Source hashes independently rechecked
- [ ] Database reopened after clean shutdown
- [ ] Backup created
- [ ] Backup restored into a separate location
- [ ] Known-answer checks pass against restored copy
- [ ] Corrupted/truncated test fixture is rejected rather than silently accepted

## Deviations

Document every deviation from the validated procedure and its disposition.

## Decision

- [ ] PASS — approved for the stated research workflow and dataset version
- [ ] FAIL — not approved

Scope/limitations:

Validator signature/name:
Date:
