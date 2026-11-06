# DNA Synthesizer Output Capsule Specification v1.0

Status: Draft.

## 1. Scope

Defines software identity and job binding for the disposable output capsule. It does not define capsule chemistry, materials, dimensions, sealing, sterility, pressure ratings or mechanical interfaces.

## 2. Cardinality

- One synthesis job MUST reserve exactly one output capsule.
- One output capsule MUST be assigned to no more than one job.
- A capsule assigned to a job cannot be reassigned to another job in v1.0, including after failed or aborted jobs.
- A new job requires a new available capsule.
- Output Capsule and Reagent Cartridge are separate components with separate IDs, inventory and lifecycle.

## 3. Capsule identity

| Field | Required | Description |
|---|---|---|
| capsule_id | Yes | Globally unique identifier |
| capsule_type | Yes | Fixed value output |
| device_id | Yes | Assigned device |
| job_id | Once reserved | Sole job binding |
| status | Yes | Capsule lifecycle state |
| assigned_at | Once reserved | UTC timestamp |
| sealed_at | When sealed | UTC timestamp |
| result_status | After job | Result status reference |

## 4. Suggested states

AVAILABLE → RESERVED → PRESENT → IN_USE → SEALED → COLLECTED

Exception states: REJECTED, CONTAMINATION_SUSPECTED, DAMAGED, UNKNOWN.

Physical presence and identity must be verified by a validated device mechanism before execution. Software reservation alone does not prove the correct capsule is physically installed.

## 5. Required checks

Before execution, verify capsule eligibility, exclusive assignment, output type, device/job/package binding and physical identity/presence.

After execution, record terminal job status and sealing/collection state only when supported by device evidence. Preserve capsule-to-job traceability. Capsule status alone does not establish product quality or sterility.

## 6. Reagent cartridge separation

Reagent cartridges use a separate entity and contract. Do not reuse capsule_id for a reagent cartridge or infer reagent identity, compatibility or remaining capacity from output capsule metadata.
