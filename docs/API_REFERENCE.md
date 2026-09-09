# API Reference Guide

This document provides a comprehensive reference for the public Rust APIs exposed by `scies-bio-db`.

---

## 1. Top-Level Entry Point: `BioDb`

`BioDb` is the central coordination handle for all storage, metadata, and ingestion operations.

```rust
use scies_bio_db::BioDb;
use std::path::Path;
```

### Methods

#### `pub fn open(db_dir: &Path) -> Result<Self>`
Opens an existing database directory, or initializes a new one with full SQLite schemas and sequence directories.
- **Errors**: Returns `BioDbError::Io` if directories cannot be created, or `BioDbError::Sql` if SQLite schema initialization fails.

#### `pub fn register_organism(&self, org: &Organism) -> Result<()>`
Inserts or updates an organism in the metadata catalog.
- **Parameters**: `org`: Instance of `scies_bio_db::model::taxonomy::Organism`.

#### `pub fn register_reference_genome(&self, genome: &ReferenceGenome) -> Result<()>`
Registers a reference genome assembly (e.g. `GRCh38.p14`, `ASM584v2`).

#### `pub fn store_sequence(&self, tax_id: u32, assembly: &str, chr: &str, sequence: &[u8]) -> Result<()>`
Encodes raw ASCII nucleotide bytes (`A`, `C`, `G`, `T`) into the 2-bit `.seq` format and writes to:
`<db_root>/seq/<tax_id>/<assembly>/<chr>.seq`

#### `pub fn fetch_sequence(&self, tax_id: u32, assembly: &str, chr: &str, start: u64, end: u64) -> Result<Vec<u8>>`
Directly reads and decodes a slice `[start, end)` from the memory-mapped sequence file.

---

## 2. Query Subsystems

### A. Genome Queries: `db.genome()`

Accessed via `let query = db.genome();`.

#### `pub fn fetch_subsequence(&self, tax_id: u32, assembly: &str, chr: &str, start: u64, end: u64) -> Result<Dna>`
Fetches sequence slice and wraps it into a verified `scies_bio_th::Dna` type.

#### `pub fn fetch_full_sequence(&self, tax_id: u32, assembly: &str, chr: &str) -> Result<Dna>`
Fetches the entire chromosome sequence from index 0 to $N$.

#### `pub fn list_chromosomes(&self, assembly: &str) -> Result<Vec<Chromosome>>`
Queries the SQLite catalog for all registered chromosomes belonging to the assembly.

---

### B. Gene & Annotation Queries: `db.genes()`

Accessed via `let query = db.genes();`.

#### `pub fn find_by_symbol(&self, symbol: &str) -> Result<Gene>`
Searches the `genes` table by official gene symbol across all registered organisms (e.g. `"INS"`, `"TP53"`).

#### `pub fn find_by_symbol_scoped(&self, tax_id: u32, symbol: &str) -> Result<Gene>`
Finds a gene by taxonomy ID and symbol, strictly isolating species lookups (e.g. `9606` for Human, `10090` for Mouse).

#### `pub fn find_by_id(&self, gene_id: &str) -> Result<Gene>`
Searches by Ensembl or HGNC gene identifier (e.g. `"ENSG00000141510"`).

#### `pub fn find_by_id_scoped(&self, gene_id: &str, assembly: &str) -> Result<Gene>`
Searches by gene identifier strictly scoped to a target reference genome assembly.

#### `pub fn find_by_region(&self, tax_id: u32, assembly: &str, chr: &str, start: u64, end: u64) -> Result<Vec<Gene>>`
Finds all genes whose genomic span overlaps the 0-based half-open range `[start, end)` scoped to `(tax_id, assembly)`.

#### `pub fn transcripts_of(&self, gene_id: &str) -> Result<Vec<Transcript>>`
Finds all transcripts associated with a given gene identifier.

---

### C. Protein Queries: `db.proteins()`

Accessed via `let query = db.proteins();`.

#### `pub fn find_by_accession(&self, accession: &str) -> Result<ProteinEntry>`
Retrieves a protein record by primary UniProtKB accession (e.g. `"P04637"`).

#### `pub fn find_by_gene(&self, gene_symbol: &str) -> Result<Vec<ProteinEntry>>`
Retrieves all protein records matching a given gene symbol.

#### `pub fn domains_of(&self, accession: &str) -> Result<Vec<Domain>>`
Retrieves protein domain annotations (Pfam, InterPro) for a protein entry.

---

### D. Variant Queries: `db.variants()`

Accessed via `let query = db.variants();`.

#### `pub fn find_by_id(&self, variant_id: &str) -> Result<VariantRecord>`
Retrieves a variant record by dbSNP rsID (e.g. `"rs1042522"`).

#### `pub fn find_by_region(&self, tax_id: u32, assembly: &str, chr: &str, start: u64, end: u64) -> Result<Vec<VariantRecord>>`
Retrieves all variants located within coordinate range `[start, end)` scoped to `(tax_id, assembly)`.

#### `pub fn find_by_gene(&self, gene_symbol: &str) -> Result<Vec<ClinicalVariant>>`
Retrieves ClinVar clinical interpretations for a given gene.

---

## 3. Ingestion Pipelines

All ingestion pipelines are atomic and require explicit taxonomy and assembly scoping:

### A. FASTA Ingestion (`FastaIngestor`)
Streams raw FASTA files, losslessly encodes 2-bit `.seq v2` files, and updates chromosome tables:
```rust
let chromosomes_ingested = db.ingest_fasta(Path::new("hg38.fa"), 9606, "GRCh38")?;
```

### B. GFF3 Annotation Ingestion (`Gff3Ingestor`)
Parses GFF3 features and populates genes, transcripts, and exon tables with validated parent-child relationships:
```rust
let features_ingested = db.ingest_gff3(Path::new("annotations.gff3"), 9606, "GRCh38")?;
```

### C. VCF Variant Ingestion (`VcfIngestor`)
Streams VCF files and loads variant records into the catalog:
```rust
let variants_ingested = db.ingest_vcf(Path::new("variants.vcf"), 9606, "GRCh38")?;
```

---

## 4. In-Memory Accelerators & Indexes

### `IntervalIndex`
Genomic interval binary-search index for fast spatial feature lookups:
```rust
use scies_bio_db::index::{IntervalIndex, IndexEntry};

let mut idx = IntervalIndex::new();
idx.insert(IndexEntry {
    chr: "chr11".to_string(),
    start: 2159778,
    end: 2161209,
    payload_id: "ENSG00000254647".to_string(),
});
idx.build();

let hits = idx.find_overlapping("chr11", 2160000, 2160500);
```

### `KmerIndex`
Rolling hash k-mer index for instant sequence pattern matching and sequence search.

---

## 5. Error Handling (`BioDbError`)

Every method returns a strongly-typed `scies_bio_db::Result<T>`:

```rust
pub enum BioDbError {
    Io(std::io::Error),
    Sql(rusqlite::Error),
    Bio(scies_bio_th::BioError),
    Encoding(String),
    NotFound(String),
    InvalidCoordinate { chr: String, start: u64, end: u64 },
    CorruptData(String),
    ChecksumMismatch { expected: u32, found: u32 },
    UnsupportedVersion(u16),
    PathTraversal(String),
}
```
