# Synthesis Integration Module (SIM) — Architecture v1.0

Status: Draft. Parent system: SCIES-BIO-DB / GENESIS ecosystem.

## 1. Purpose

SIM provides a controlled software boundary between genome-reference retrieval and DNA Synthesizer job intake. It converts an explicitly selected reference interval into an immutable, provenance-bearing execution package. It does not design or edit sequences and does not implement synthesis chemistry.

## 2. Principles

1. Reference-only source: resolve genome records from SCIES-BIO-DB; no user genome-file ingestion in v1.0.
2. No sequence editing: extraction and optional reverse-complement orientation are the only permitted sequence transformations.
3. Immutable provenance: preserve taxon, assembly, contig, interval, strand, source identity and integrity evidence.
4. Fail closed: unresolved references, invalid intervals, ambiguous bases, integrity mismatches, missing authorization or capsule conflicts block readiness.
5. Separation of concerns: genome storage, job orchestration, device communications and physical synthesis remain separate.
6. Local start confirmation: package receipt never starts the instrument.
7. Auditable lifecycle: transitions and acknowledgements are timestamped and correlated by IDs.

## 3. Components

### Genome Reference Adapter
Calls the existing SCIES-BIO-DB API; resolves tax_id, assembly, contig/sequence ID, contig length and reference identity; retrieves the requested region; invokes available integrity verification. It MUST use supported database APIs rather than query SQLite internals directly.

### Coordinate and Extraction Service
Normalizes source coordinates to 0-based, half-open coordinates; validates interval bounds and length; applies reverse-complement orientation only when explicitly requested; records original and normalized coordinates.

### Sequence Contract Validator
Enforces uppercase A/C/G/T for v1.0 automatic intake; rejects ambiguous IUPAC symbols, U, gaps and malformed content; computes canonical sequence SHA-256. It performs no substitutions, trimming, codon changes, repair or other sequence editing.

### Job and Authorization Service
Creates immutable job and sequence-version records; tracks policy/screening and technical-review references; prevents package readiness until required review gates pass; keeps authorization separate from provenance.

### Package and Integrity Service
Produces schema-valid packages, deterministic manifest integrity metadata and (in production) digital signatures backed by managed keys. Rejects unsupported contract versions.

### Capsule Registry
Reserves one output capsule for one job; prevents duplicate/concurrent assignment; tracks capsule identity separately from reagent cartridges.

### Device Gateway
Transfers packages using a transport adapter selected by the product team; receives acknowledgements and events; implements idempotency and replay protection; requires device-local confirmation. It does not expose database credentials to the device.

### Result and Audit Service
Stores device/job/package/capsule IDs, timestamps, transitions, verification evidence and terminal status. It distinguishes completed, failed, aborted, interrupted and unknown.

## 4. Proposed module layout

    synthesis-integration/
    ├── src/
    │   ├── genome_source/       # SCIES-BIO-DB adapter
    │   ├── coordinates/         # normalization and bounds checks
    │   ├── extraction/          # reference region extraction
    │   ├── sequence/            # validation and hashing
    │   ├── job/                 # lifecycle and authorization
    │   ├── package/             # manifest, integrity, signing
    │   ├── capsule/             # output capsule registry
    │   ├── device/              # transport and acknowledgements
    │   └── result/              # result and audit persistence
    ├── schemas/
    └── tests/

This is a proposed logical layout, not a claim that these source files already exist.

## 5. Data flow

1. GENESIS selects a genome reference and interval.
2. SIM resolves the exact assembly and contig in SCIES-BIO-DB.
3. SIM verifies the source and extracts the interval.
4. SIM normalizes coordinates and applies the explicitly selected strand orientation.
5. SIM validates alphabet, interval length and provenance.
6. Separate screening, technical review and authorization gates determine whether the job may proceed.
7. SIM reserves one output capsule and creates an immutable job.
8. SIM creates and integrity-protects the execution package.
9. Device Gateway transfers the package.
10. Device validates it and acknowledges receipt.
11. Device waits for local physical confirmation.
12. Events and final result are recorded against the same job/package/capsule IDs.

## 6. Logical interfaces

- resolve_reference(source_id, assembly, seq_id)
- extract_region(reference, start, end, strand)
- validate_sequence(sequence, provenance)
- create_job(sequence_ref, authorization_ref, capsule_id, idempotency_key)
- build_package(job_id)
- submit_package(device_id, package_id)
- record_device_event(device_id, job_id, event)
- record_result(job_id, result)

Transport and exact API framework remain undecided.

## 7. Security and reliability

Authenticate GENESIS, SIM and device identities. Authorize operations using least privilege. Protect data in transit and at rest. Sign production packages and maintain key rotation/revocation. Enforce idempotency and replay protection. Use append-only or tamper-evident audit events. Do not place database credentials or signing private keys on the device. Treat malformed, partial, stale or unverifiable packages as blocked. Define recovery behavior for power and communication loss before hardware validation.

## 8. SCIES-BIO-DB integration

Existing genome/reference APIs and SequenceRef-like metadata are the source adapter boundary. SIM consumes stable public APIs rather than reading SQLite tables or storage files directly. The compact .seq representation is a storage concern; SIM decodes through the supported DB layer and creates canonical nucleotide text before contract hashing.

The existing 2-bit encoding is not a device command format. Ambiguous-base sidecar data must not be silently coerced into A/C/G/T.

Do not insert synthesis-job tables into the genome core until persistence ownership and API boundaries are agreed. Initial job metadata may live in a separate module/service.

## 9. Open decisions

- Exact SCIES-BIO-DB API method and coordinate convention for each retrieval path.
- Job persistence technology and ownership.
- Authentication/signature format and key lifecycle.
- Device transport protocol.
- Screening service interface and approval authority.
- Supported oligo length and chemistry compatibility.
- Physical capsule identification and verification.
- Result verification evidence and telemetry contract.
