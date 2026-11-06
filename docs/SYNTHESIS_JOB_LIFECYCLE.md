# DNA Synthesizer Job Lifecycle v1.0

Status: Draft.

## 1. States

| State | Meaning | Allowed next states |
|---|---|---|
| CREATED | Job request recorded | REFERENCE_RESOLVED, REJECTED |
| REFERENCE_RESOLVED | Exact assembly and contig resolved | EXTRACTED, REJECTED |
| EXTRACTED | Region extracted from reference | VALIDATED, REJECTED |
| VALIDATED | Coordinates, alphabet, length and hashes pass | REVIEW_PENDING, REJECTED |
| REVIEW_PENDING | Awaiting policy and technical review | APPROVED, REJECTED |
| APPROVED | Required authorization recorded | PACKAGE_READY, REVOKED |
| PACKAGE_READY | Immutable package and capsule binding created | DEVICE_RECEIVED, CANCELLED |
| DEVICE_RECEIVED | Device received package | DEVICE_VALIDATED, REJECTED |
| DEVICE_VALIDATED | Device verified package and binding | AWAITING_LOCAL_CONFIRMATION, REJECTED |
| AWAITING_LOCAL_CONFIRMATION | Waiting for local user confirmation | EXECUTING, CANCELLED |
| EXECUTING | Device reports execution started | VERIFYING, FAILED, ABORTED, INTERRUPTED |
| VERIFYING | Result evidence being evaluated | COMPLETED, FAILED, UNKNOWN |
| COMPLETED | Required terminal evidence persisted | Terminal |
| FAILED | Execution or verification failed | Terminal; new job required |
| ABORTED | Explicitly aborted | Terminal |
| INTERRUPTED | Execution interrupted; outcome not established | VERIFYING, UNKNOWN |
| UNKNOWN | Outcome cannot be established | Terminal pending controlled reconciliation |
| REJECTED | Validation or authorization gate failed | Terminal |
| REVOKED | Authorization withdrawn before execution | Terminal |
| CANCELLED | Cancelled before execution | Terminal |

Implementations may add states but MUST NOT weaken gates or reinterpret terminal failures as success.

## 2. Required gates

### Reference gate
Exact reference, assembly and contig must resolve. Interval must be in bounds and coordinate semantics explicit.

### Sequence gate
Sequence must be canonical uppercase A/C/G/T, non-empty, and length must equal normalized interval length. Source and sequence integrity checks must pass. No editing is permitted.

### Review and authorization gate
A format check is not authorization to synthesize. Required biological safety/policy screening, technical feasibility review and human approval must have explicit references and outcomes. Missing or indeterminate review blocks package creation.

### Capsule gate
One eligible output capsule must be reserved exclusively for the job. Capsule ID and job ID must be bound before package release.

### Device gate
Device identity, package version, manifest integrity, authorization, capsule binding and replay/idempotency checks must pass. Device waits for local confirmation.

### Result gate
Terminal success requires device-reported completion and defined result-verification evidence. Missing or contradictory evidence results in UNKNOWN or FAILED, never COMPLETED.

## 3. Event record

Every transition SHOULD record event_id, job_id, package_id when available, previous_state, new_state, actor_type, actor_id, device_id when applicable, UTC RFC 3339 timestamp, reason_code and evidence references.

Events should be append-only or tamper-evident.

## 4. Retry and idempotency

- Repeating a create request with the same idempotency key and same body returns the original job.
- Reusing an idempotency key with a different body is a conflict.
- Device receipt and acknowledgement retries must not create a second execution.
- A job that entered EXECUTING must not be restarted through blind retry.
- After uncertain interruption, reconcile device state and evidence before further action. If safe state cannot be established, retain UNKNOWN and require a new controlled workflow.

## 5. Cancellation and revocation

Cancellation or authorization revocation before execution prevents the device from starting. After execution begins, a cancellation request is not proof that execution stopped. Final state must be based on device evidence.
