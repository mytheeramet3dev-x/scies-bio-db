# scies-bio-db v0.2 Architecture & Phase DB-1 Design Specification

## 1. Overview and Mission

`scies-bio-db` is the genomic and biological annotation engine of the SCIES ecosystem (`scies-bio-th`, `scies-bio-db`, `Humanz`). Its mission is to deliver deterministic, memory-efficient storage and sub-millisecond retrieval of reference chromosomes, gene annotations, and longevity intervention loci, without depending on external database servers or non-reproducible web queries.

Phase DB-1 transitions `scies-bio-db` from a theoretical blueprint into a hardened, production-ready storage engine for real human genomic data (T2T-CHM13v2.0 and GRCh38).

---

## 2. Architecture and Storage Layout

`scies-bio-db` uses a hybrid architecture separating high-volume raw sequence payloads from structured metadata:

```text
scies_bio_db::
├── api            — BioDb public facade (single point of entry)
├── model          — Taxonomy, Assembly, SequenceRef, GeneFeature, LongevityTarget, LongevityLocus
├── storage        — 2-bit sequence files (.seq v2) + SQLite catalog (meta.db)
├── index          — Coordinate, K-mer, and Interval lookup indices
├── ingest         — FASTA, GFF3/GTF, and VCF streaming importers
├── query          — Typed query handles (GenomeQuery, GeneQuery, etc.)
└── error          — BioDbError (thiserror-based typed error hierarchy)
```

### On-Disk Directory Layout

```text
<db_dir>/
├── meta.db                         # Embedded SQLite catalog (WAL mode, busy timeout 5s)
└── seq/
    └── <tax_id>/                   # NCBI Taxonomy ID (e.g. 9606 for Homo sapiens)
        └── <assembly_name>/        # Reference assembly (e.g. "T2T-CHM13v2.0")
            ├── chr1.seq            # Memory-mapped .seq v2 file
            ├── chr2.seq
            └── ...
```

---

## 3. Storage Engine: `.seq v2` Format

### 3.1 Binary Specification
A `.seq v2` file comprises:
1. **64-byte Header**:
   - `[0..8]`: Magic bytes `b"SCIESSEQ"`
   - `[8..10]`: Version `u16` in little-endian (`2`)
   - `[10..12]`: Reserved
   - `[12..20]`: Total sequence length `u64`
   - `[20..28]`: Packed 2-bit payload byte length `u64`
   - `[28..36]`: Ambiguity record count `u64`
   - `[36..40]`: IEEE 802.3 CRC32 checksum `u32` of packed payload and sidecar
   - `[40..64]`: Zero padding for future extension
2. **2-bit Packed Sequence Payload**:
   - `A = 00`, `C = 01`, `G = 10`, `T = 11` (packed 4 bases per byte, MSB first)
3. **Lossless Ambiguity Sidecar**:
   - Sorted list of `(u64 position, u8 ASCII code)` records representing non-ACGT IUPAC bases (`N`, `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, and Uracil `U`).

### 3.2 Integrity and Safety Invariants
- **Atomic Persistence**: Files are written to `.tmp.<random>` siblings, flushed, synced, and renamed atomically to prevent partial writes.
- **Path Isolation**: All taxonomy, assembly, and chromosome string parameters pass through `sanitize_component` to prevent directory traversal attacks (`..`, `/`, `\`).
- **Memory-Mapped Verification**: Cached memory maps validate magic, version, lengths, and CRC32 before servicing slice reads.

---

## 4. SQLite Metadata Store & Provenance

### 4.1 Core Schema

```sql
-- Reference Assemblies
CREATE TABLE IF NOT EXISTS assemblies (
    tax_id        INTEGER NOT NULL,
    assembly_name TEXT    NOT NULL,
    source        TEXT    NOT NULL,
    created_at    TEXT    NOT NULL,
    notes         TEXT,
    PRIMARY KEY (tax_id, assembly_name)
);

-- Reference Chromosome & Contig Sequences
CREATE TABLE IF NOT EXISTS sequences (
    tax_id        INTEGER NOT NULL,
    assembly_name TEXT    NOT NULL,
    seq_id        TEXT    NOT NULL,
    length_bp     INTEGER NOT NULL,
    path          TEXT    NOT NULL,
    sha256        TEXT    NOT NULL,
    source_uri    TEXT,
    PRIMARY KEY (tax_id, assembly_name, seq_id)
);

-- Gene Annotations
CREATE TABLE IF NOT EXISTS genes (
    gene_id       TEXT    NOT NULL,
    assembly_name TEXT    NOT NULL,
    symbol        TEXT    NOT NULL,
    name          TEXT    NOT NULL,
    biotype       TEXT    NOT NULL,
    chr           TEXT    NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        TEXT    NOT NULL,
    tax_id        INTEGER NOT NULL,
    feature_type  TEXT    NOT NULL DEFAULT 'gene',
    source        TEXT,
    PRIMARY KEY (gene_id, assembly_name)
);

-- Curated Longevity Intervention Targets
CREATE TABLE IF NOT EXISTS longevity_targets (
    symbol           TEXT    NOT NULL,
    tax_id           INTEGER NOT NULL,
    priority         INTEGER NOT NULL,
    category         TEXT    NOT NULL,
    default_assembly TEXT,
    notes            TEXT,
    PRIMARY KEY (tax_id, symbol)
);
```

### 4.2 Cryptographic Provenance
Every reference sequence ingested computes its SHA-256 hash using streaming computation:
```rust
let sha256 = scies_bio_th::seq::compute_sha256_hex(&sequence);
```
Stored in `sequences.sha256`, this hash is verified on demand via `BioDb::verify_sequence_integrity(tax_id, assembly, seq_id)` to catch bit-rot or storage tampering.

---

## 5. Ingestion Pipeline

### 5.1 FASTA Importer
- Streams records with memory-bounded buffers.
- Normalizes sequence IDs by stripping whitespace and extra header tokens.
- Encodes into `.seq v2` and registers `Assembly` and `SequenceRef` in SQLite within transactions.

### 5.2 GFF3 / GTF Importer
- Resolves gene symbols across standard attributes: `gene_name`, `Name`, `gene`, and `symbol`.
- Validates coordinates: invalid or inverted coordinates (`start > end`) are tracked in `report.skipped_count` rather than crashing the ingestion batch.
- Automatically normalizes coordinates into 0-based half-open intervals `[start, end)`.

---

## 6. Curated Longevity Catalog

`scies-bio-db` seeds 10 essential human longevity intervention targets on initialization:

| Symbol | Category | Priority | Default Assembly | Mechanism / Role |
| :--- | :--- | :---: | :--- | :--- |
| **TERT** | `telomere` | 1 | T2T-CHM13v2.0 | Catalytic reverse transcriptase for telomere elongation |
| **TERC** | `telomere` | 1 | T2T-CHM13v2.0 | Essential template RNA subunit of telomerase |
| **SIRT6** | `dna_repair` | 1 | T2T-CHM13v2.0 | NAD+-dependent deacetylase; chromatin & DSB repair |
| **TP53** | `tumor_suppressor` | 1 | T2T-CHM13v2.0 | Safeguards reprogramming and suppresses aberrant clones |
| **ATM** | `dna_repair` | 2 | T2T-CHM13v2.0 | Master orchestrator of double-strand break response |
| **CDKN2A**| `tumor_suppressor` | 2 | T2T-CHM13v2.0 | p16INK4a / p14ARF; cellular senescence driver |
| **FOXO3** | `metabolic` | 2 | T2T-CHM13v2.0 | Longevity-associated transcription factor (GWAS) |
| **KLOTHO**| `metabolic` | 2 | T2T-CHM13v2.0 | Anti-aging circulating endocrine factor |
| **LMNA** | `proteostasis` | 2 | T2T-CHM13v2.0 | Nuclear envelope integrity; mutated in progeria |
| **PARP1** | `dna_repair` | 2 | T2T-CHM13v2.0 | Base excision repair and early DNA damage signaling |

### 6.1 `LongevityLocus` Model
Combines:
- `gene: GeneFeature`
- `target: Option<LongevityTarget>`
- `sequence: scies_bio_th::Dna`
- `window_start: u64` / `window_end: u64`
- `provenance_sha256: String`
- `tax_id: u32` / `assembly: String`

---

## 7. Downstream Integration

`LongevityLocus` feeds directly into `scies-bio-th` algorithms:
1. **CRISPR Guide Design**: `scies_bio_th::aging::design_real_guides(&locus.sequence, strand, locus.window_start)`
2. **Promoter & CpG Island Scanning**: `scies_bio_th::epigenome::methylation::find_cpg_islands(&locus.sequence, ...)`
3. **Genome Architecture Auditing in `Humanz`**: Providing verified real genomic loci for whole-genome simulation.

---

## 8. Definition of Done Checklist (Phase DB-1)

- [x] `cargo test --release`: 26 passed, 0 failed (100% pass rate)
- [x] `cargo check --release --all-targets`: 0 errors, 0 warnings
- [x] `BioDb::open(...)` with FASTA and GFF3 ingestion working
- [x] Typed query layer: `get_sequence`, `get_gene_by_symbol`, `list_longevity_targets`, `get_longevity_locus`
- [x] Cryptographic SHA-256 provenance on all sequences and assemblies
- [x] Runnable example: `examples/ingest_and_query_longevity.rs`
- [x] README and architecture documentation updated
- [x] Downstream compatibility with `scies-bio-th` and `Humanz` preserved
