# DNA Synthesizer Data Contract v1.0

Status: Draft specification for software integration and digital prototyping.

## 1. Locked decisions

- Genome source: SCIES-BIO-DB reference data only.
- Sequence editing: none. GENESIS may resolve a reference, select a region and extract it; it MUST NOT alter nucleotide content.
- User-imported genome files and manually edited sequences are out of scope for v1.0.
- The device receives an immutable, versioned genome-derived sequence package.
- The device validates the package and waits for explicit local confirmation before execution.
- Exactly one output capsule is assigned to one synthesis job. Output capsules and reagent cartridges are separate entities.

## 2. System boundary

SCIES-BIO-DB is the authoritative genome-reference source. GENESIS Computer is the orchestration client. The Synthesis Integration Module (SIM) resolves reference metadata, extracts the requested region, validates provenance and creates an immutable job package. The device does not query SCIES-BIO-DB directly.

SCIES-BIO-DB → GENESIS Computer / SIM → Execution Package → DNA Synthesizer → Result Record

The package is a data and authorization envelope, not a chemical synthesis recipe or validated instrument process program.

## 3. Canonical genome coordinates

SIM MUST normalize intervals to 0-based, half-open coordinates [start, end). It MUST retain the source coordinate convention and original coordinates when the upstream source uses another convention.

Validation requirements:
- start >= 0
- end > start
- end - start equals extracted sequence length
- interval is within the referenced contig length
- assembly and contig identifiers resolve unambiguously

Strand is explicit: + or -. For -, SIM may derive the reverse complement as an orientation operation. It MUST record this operation and MUST NOT otherwise edit bases.

## 4. Sequence rules

The v1.0 canonical alphabet is uppercase DNA A, C, G and T. Ambiguous IUPAC symbols, U, whitespace inside the sequence, gaps and other characters MUST be rejected from automatic job creation. They MUST NOT be silently converted, removed or guessed.

Sequence content is immutable after extraction. If the selected reference, interval, strand or extraction result changes, create a new sequence/job identity and hash. Never overwrite an existing job sequence.

Passing this data contract is not evidence that a sequence is biologically safe, chemically compatible, manufacturable or approved. Separate policy screening, technical feasibility and human approval gates are mandatory before execution.

## 5. Integrity and identity

- SHA-256 is used for canonical sequence and manifest integrity.
- Sequence hash input is UTF-8 bytes of the uppercase canonical sequence, without whitespace or trailing newline.
- Manifest hash is calculated over deterministic canonical serialization with the manifest hash field omitted. The canonical JSON algorithm and version MUST be documented.
- A hash detects content changes; it does not authenticate the sender. Production deployments SHOULD sign manifests and verify signer trust independently.
- Job creation and device submission MUST be idempotent.

## 6. Execution package

The package MUST include contract version, package/job IDs, immutable sequence identity and content, genome provenance, output capsule identity, authorization reference and integrity metadata. It MUST NOT include unvalidated chemistry parameters or instrument-level synthesis instructions.

The device MUST verify schema, package integrity, authorization validity, supported contract version, sequence constraints, capsule/job binding and replay/idempotency state. After successful receipt and validation, it MUST enter AWAITING_LOCAL_CONFIRMATION. Receipt alone MUST NOT start the instrument.

## 7. Job and capsule rules

One job binds to exactly one output capsule. A capsule MUST NOT be assigned to multiple jobs. A job MUST NOT be reported completed unless the device has returned a terminal success result and the result record has been persisted. Failure, abort, interruption and unknown state are distinct outcomes.

## 8. Versioning

Consumers MUST reject unsupported major contract versions. Additive optional fields may be introduced in a compatible minor revision. Changes to coordinate semantics, sequence canonicalization, hash inputs or required identity fields require a major version change.

## 9. Non-goals

This contract does not specify synthesis chemistry, reagent formulation, fluidics, temperature, cycle timing, yield, purity, oligo length limits, device mechanics, capsule material, sterilization or a validated manufacturing protocol. These require separate engineering and validation documents.

## 10. Related documents

- SYNTHESIS_INTEGRATION_ARCHITECTURE.md
- SYNTHESIS_JOB_LIFECYCLE.md
- SYNTHESIS_CAPSULE_SPEC.md
- JSON schemas in the schemas directory.
