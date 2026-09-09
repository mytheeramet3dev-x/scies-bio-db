# Storage Engine & System Architecture

`scies-bio-db` is an enterprise-grade, domain-specific database engine optimized for biological sequences and genomic metadata. It addresses the fundamental asymmetry in biological data:

1. **Sequence Payload**: Massive, append-mostly, dense, linear nucleotide arrays (billions of bases) requiring contiguous memory mapping and fast random slice access.
2. **Annotation & Metadata**: Highly structured, relational, searchable biological records (organisms, genes, exons, transcripts, variants, proteins, phenotypes) requiring ACID compliance and multi-parameter indexed queries.

To achieve maximum performance with minimal disk footprint, `scies-bio-db` uses a **hybrid dual-engine architecture**:

```text
                                  ┌────────────────────────┐
                                  │   BioDb Public API     │
                                  │ (Query & Ingestion Hub)│
                                  └───────────┬────────────┘
                                              │
                    ┌─────────────────────────┴─────────────────────────┐
                    ▼                                                   ▼
       ┌────────────────────────┐                          ┌────────────────────────┐
       │     SeqStore (I/O)     │                          │   SqliteStore (Meta)   │
       │  2-Bit Mmap Sequences  │                          │ Embedded SQLite Engine │
       └────────────┬───────────┘                          └────────────┬───────────┘
                    │                                                   │
     ┌──────────────┴──────────────┐                     ┌──────────────┴──────────────┐
     ▼                             ▼                     ▼                             ▼
[.seq Files]                [Zero-Copy Mmap]       [meta.db Tables]            [Relational Indices]
- 4 bases / byte            - Direct virtual       - organisms                 - gene symbol index
- 8-byte LE header            address translation  - chromosomes               - interval ranges
- Partitioned by TaxID      - OS page caching      - genes & transcripts       - UniProt accession
```

---

## 1. Storage Layers

### A. The 2-Bit Sequence Engine (`SeqStore`)
- **Format**: Custom `.seq` flat binary files (see [FORMAT_SPEC.md](FORMAT_SPEC.md) for full binary layout).
- **Compression**: Standard ASCII nucleotide representations require 1 byte per character (`0x41` for 'A', etc.). `SeqStore` packs 4 bases into 1 byte (2 bits per base), reducing disk footprint by exactly **75%**.
- **I/O Subsystem**: Reads utilize `memmap2::MmapOptions` to establish zero-copy memory maps directly into the kernel page cache. Random slices (`[start..end)`) are read with simple pointer arithmetic, bypassing user-space buffer allocations.

### B. The Relational Metadata Engine (`SqliteStore`)
- **Database Engine**: Embedded SQLite 3 via `rusqlite` (compiled with the `bundled` feature for zero C/system dependency issues).
- **Concurrency**: The underlying connection is guarded by an `Arc<Mutex<Connection>>` for safe multi-threaded query execution across threads.
- **Relational Tables**:
  - `organisms`: Taxonomic lineage, common names, NCBI Taxonomy IDs.
  - `reference_genomes`: Assembly names, patch status, release years.
  - `chromosomes`: Chromosome names, sequence lengths, karyotype classifications.
  - `genes`: Ensembl/HGNC IDs, gene symbols, biotypes, genomic coordinates.
  - `transcripts`: Transcript IDs, biotypes, exon counts.
  - `exons`: Genomic boundaries, exon numbers, phases.
  - `protein_entries`: UniProt accessions, subcellular locations, review status.
  - `variant_records`: dbSNP identifiers, alleles, clinical interpretations.

---

## 2. Directory Hierarchy & Taxonomy Partitioning

Files are organized deterministically using NCBI Taxonomy IDs and Assembly releases:

```text
<db_root>/
├── meta.db                         # SQLite catalog with indexed schemas
└── seq/                            # Sequence payload root
    ├── 9606/                       # NCBI TaxID: Homo sapiens
    │   └── GRCh38.p14/             # Assembly name
    │       ├── chr1.seq            # Chromosome sequence files
    │       ├── chr2.seq
    │       └── chr11_INS.seq
    ├── 562/                        # NCBI TaxID: Escherichia coli
    │   └── ASM584v2/
    │       └── lacZ.seq
    └── 2562234/                    # NCBI TaxID: Panthera spelaea (Cave Lion)
        └── PanSpe_Yakutia_v1.0/
            └── ExUtero_Clone_Cub_01.seq
```

### Why Taxonomy-First Partitioning?
1. **Name Collision Elimination**: Different species share chromosome designations (`chr1`, `chrX`, `MT`). Partitioning by numerical Taxonomy ID isolates species namespace completely.
2. **Multi-Assembly Coexistence**: Multiple reference assemblies (e.g. `GRCh37` and `GRCh38`) reside cleanly under the same TaxID.
3. **OS Filesystem Optimization**: Prevents directory bloat when storing tens of thousands of chromosomes or viral segments.

---

## 3. Concurrency & Memory Model

- **Reader Threads**: Multiple concurrent readers can invoke `fetch_sequence` simultaneously. Because each `.seq` file is memory-mapped in read-only mode, the operating system kernel handles concurrent page reads without lock contention.
- **Writer Operations**: Writing or updating sequences performs append-only writes via `BufWriter<File>` to temporary files before atomic file renaming, preventing torn reads.
- **Thread Safety**: All query wrappers (`GenomeQuery`, `GeneQuery`, `ProteinQuery`, `VariantQuery`) hold cheap `Arc` references to `SeqStore` and `SqliteStore`, making `BioDb` safe to clone and share across worker threads (`Send + Sync`).
