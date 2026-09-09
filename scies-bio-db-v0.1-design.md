# scies-bio-db — v0.1 Design Document

> A biological database engine for all living organisms, starting from *Homo sapiens*.
> ส่วนหนึ่งของซีรีส์ `scies` ทำงานคู่กับ `scies-bio-th` โดยตรง

---

## 1. วิสัยทัศน์และเป้าหมาย

`scies-bio-db` คือ **ฐานข้อมูลชีววิทยาครอบจักรวาล** ที่เก็บข้อมูลของสิ่งมีชีวิตทุกสปีชีส์บนโลก
โดยเริ่มต้นจาก **มนุษย์ (*Homo sapiens*, Taxonomy ID: 9606)** เป็นแกนหลักก่อน
จากนั้นขยายไปยัง Model Organisms และสิ่งมีชีวิตอื่น ๆ ตามลำดับ

### หลักการออกแบบ

- **Separation of Concerns**: แยก Computation (`scies-bio-th`) ออกจาก Storage (`scies-bio-db`) ชัดเจน
- **Type Safety**: ใช้ Rust type system ป้องกัน coordinate confusion (0-based vs 1-based, strand direction)
- **Fast Access**: Memory-mapped sequence files สำหรับ Random access ลำดับเบสโดยไม่โหลดทั้งโครโมโซมเข้า RAM
- **Layered Data**: ข้อมูลแบ่งเป็น 6 ชั้น ตั้งแต่ Taxonomy ไปจนถึง Systems Biology
- **scies-bio-th First**: ใช้ Types และ Parsers ของ `scies-bio-th` เป็น Data Model กลางทุกที่

---

## 2. ข้อมูลที่จัดเก็บ (6 Layers)

### Layer 1: Taxonomy & Organism Registry
ทะเบียนสิ่งมีชีวิตทุกตัวในระบบ — ทุก entity ในชั้นอื่นอ้างอิงมาที่ `tax_id`

```
Organism { tax_id: u32, scientific_name, common_name, lineage: Lineage }
Lineage  { domain, kingdom, phylum, class, order, family, genus, species }
```

**มนุษย์**: `Organism { tax_id: 9606, scientific_name: "Homo sapiens", ... }`

### Layer 2: Genomic Sequences
จีโนมอ้างอิงและลำดับเบสดิบ

```
ReferenceGenome { assembly_name, tax_id, status, release_year }
Chromosome      { name, length, assembly_name, is_mitochondrial }
[storage]       chr1.seq, chr2.seq, ...  ← 2-bit encoded, memory-mapped
```

**มนุษย์**: GRCh38/hg38 (23 chromosome pairs + MT) หรือ T2T-CHM13 (telomere-to-telomere)

### Layer 3: Gene Models & Annotations
โครงสร้างยีน, Transcript, Exon

```
Gene       { gene_id (ENSG/HGNC), symbol, name, biotype, chr, start, end, strand }
Transcript { transcript_id (ENST), gene_id, biotype, exon_count, ... }
Exon       { exon_id, transcript_id, chr, start, end, strand, exon_number }
CdsRegion  { transcript_id, chr, start, end, phase }
UtrRegion  { transcript_id, chr, start, end, is_five_prime }
```

**มนุษย์**: GENCODE v45 — ~20,000 protein-coding genes, ~60,000 total genes

### Layer 4: Proteomics & Functional Annotation
โปรตีนและคำอธิบายการทำงาน

```
ProteinEntry     { uniprot_accession, entry_name, gene_symbol, sequence, length, ... }
ProteinIsoform   { isoform_id, uniprot_accession, sequence }
Domain           { name, db (Pfam/InterPro), start_aa, end_aa }
Ptm              { position, ptm_type (Phosphorylation/Glycosylation/...) }
GoAnnotation     { go_id, term, aspect (BP/MF/CC) }
SubcellularLoc   { Nucleus | Cytoplasm | Membrane | Mitochondria | ... }
```

**มนุษย์**: UniProtKB Swiss-Prot Human Reference — ~20,400 reviewed entries

### Layer 5: Variants & Clinical Data
การกลายพันธุ์และความสัมพันธ์กับโรค

```
VariantRecord    { variant_id (rsID), chr, position, reference, alternate, type }
ClinicalVariant  { clinical_significance (P/LP/VUS/LB/B), condition, review_status }
VariantEffect    { transcript_id, effect_type (Missense/Frameshift/...), aa_change }
GwasAssociation  { snp_id, phenotype, p_value, odds_ratio }
```

**มนุษย์**: ClinVar + dbSNP + GWAS Catalog

### Layer 6: Systems & Pathways *(planned v0.4)*
วิถีชีวเคมีและเครือข่าย

```
Pathway  { pathway_id (KEGG/Reactome), name, members: Vec<gene_id> }
Reaction { reactants, products, enzymes }
Ppi      { protein_a, protein_b, score, evidence }
Grn      { tf_gene, target_gene, regulation_type }
```

---

## 3. สถาปัตยกรรม Storage

### 3.1 Sequence Storage (Binary Files)

ลำดับ DNA เก็บในรูปแบบ **2-bit encoding** เพื่อประหยัดพื้นที่ 75%:

```
A = 00,  C = 01,  G = 10,  T = 11
```

จีโนมมนุษย์ GRCh38 (~3.2 Gb raw) → **~800 MB** หลัง 2-bit encode

```
human_grch38/
├── seq/
│   ├── chr1.seq      ← 2-bit encoded, memory-mappable (mmap)
│   ├── chr2.seq
│   └── chrM.seq
└── idx/
    ├── chr1.idx      ← offset + length index
    └── ...
```

**Random access**: ดึงช่วง `chr1:1,000,000-1,001,000` โดยคำนวณ byte offset โดยตรง
ไม่ต้อง scan ทั้งไฟล์

### 3.2 Metadata Storage (SQLite)

Entity ทุกชั้น (ยกเว้น raw sequences) เก็บใน **SQLite** ผ่าน `rusqlite`:

```sql
-- ตัวอย่าง schema หลัก
CREATE TABLE organisms (tax_id INTEGER PRIMARY KEY, scientific_name TEXT, ...);
CREATE TABLE chromosomes (name TEXT, assembly_name TEXT, length INTEGER, ...);
CREATE TABLE genes (gene_id TEXT PRIMARY KEY, symbol TEXT, chr TEXT, start INTEGER, end INTEGER, strand TEXT, biotype TEXT, tax_id INTEGER);
CREATE TABLE transcripts (transcript_id TEXT PRIMARY KEY, gene_id TEXT, ...);
CREATE TABLE protein_entries (uniprot_accession TEXT PRIMARY KEY, sequence TEXT, ...);
CREATE TABLE variant_records (variant_id TEXT, chr TEXT, position INTEGER, ref TEXT, alt TEXT, ...);
CREATE TABLE clinical_variants (variant_id TEXT, clinical_significance TEXT, condition TEXT, ...);
```

### 3.3 Cache Layer (In-Memory KV Store)

KvStore (`HashMap<String, Vec<u8>>`) ใช้เป็น Cache ระหว่าง session
ในอนาคตสามารถเปลี่ยนเป็น `redb` หรือ `sled` ได้โดยไม่ต้องแก้ public API

---

## 4. Index Layer

### Genomic Interval Index
ค้นหา Gene / Feature ที่อยู่ในช่วง Coordinate ที่ระบุ

```rust
// ตัวอย่างใช้งาน
let genes = db.genes().find_by_region("chr17", 7_668_000, 7_690_000)?;
// → [Gene { symbol: "TP53", ... }]
```

Implementation: sorted `Vec<IndexEntry>` + binary search (O(log n))
อนาคต: Augmented Interval Tree สำหรับ performance ที่ดีขึ้น

### K-mer Hash Index
ค้นหาลำดับเบสแบบ exact match ด้วย rolling hash

```rust
let hits = db.genome().search_kmer(b"ATGCGATCGATCG")?;
// → [(chr: "chr17", pos: 7_123_456), ...]
```

---

## 5. Module Tree

```
scies_bio_db::
├── api/
│   └── BioDb          ← entry point สำหรับผู้ใช้ (open, genome(), genes(), ...)
├── model/
│   ├── taxonomy       ← Organism, Lineage, TaxRank
│   ├── genome         ← ReferenceGenome, Chromosome, Strand, AssemblyStatus
│   ├── gene           ← Gene, Transcript, Exon, CdsRegion, UtrRegion, Biotype
│   ├── protein        ← ProteinEntry, Domain, Ptm, GoAnnotation, SubcellularLocation
│   └── variant        ← VariantRecord, ClinicalVariant, VariantEffect
├── storage/
│   ├── seq_store      ← 2-bit encode/decode + memory-mapped file I/O
│   ├── sqlite_store   ← SQLite connection + schema init + CRUD
│   └── kv_store       ← in-memory cache
├── index/
│   ├── interval_idx   ← sorted interval search
│   └── kmer_idx       ← k-mer rolling hash lookup
├── ingest/
│   ├── fasta          ← Reference genome ingestion (wraps scies_bio_th::io::fasta)
│   ├── gff3           ← Gene annotation ingestion (wraps scies_bio_th::io::gff)
│   └── vcf            ← Variant ingestion (wraps scies_bio_th::io::vcf)
├── query/
│   ├── genome_query   ← fetch_sequence, list_chromosomes
│   ├── gene_query     ← find_by_symbol, find_by_region, transcripts_of
│   ├── protein_query  ← find_by_accession, find_by_gene, domains_of
│   └── variant_query  ← find_by_id, find_by_region, find_by_gene
└── error              ← BioDbError, Result<T>
```

---

## 6. Public API (ตัวอย่างการใช้งาน)

```rust
use scies_bio_db::BioDb;
use std::path::Path;

fn main() -> scies_bio_db::Result<()> {
    // เปิดฐานข้อมูล
    let db = BioDb::open(Path::new("./human_grch38"))?;

    // ดึงลำดับ DNA บนโครโมโซม 17 (ตำแหน่งยีน TP53)
    let dna = db.genome().fetch_sequence("chr17", 7_668_402, 7_687_550)?;
    println!("TP53 region: {} bp", dna.len());

    // ค้นหายีนด้วยชื่อ
    let tp53 = db.genes().find_by_symbol("TP53")?;
    println!("Gene: {} | Biotype: {:?}", tp53.symbol, tp53.biotype);

    // ดึง Transcript ทั้งหมดของ TP53
    let transcripts = db.genes().transcripts_of(&tp53.gene_id)?;
    println!("Transcripts: {}", transcripts.len());

    // ดึงข้อมูลโปรตีน (Canonical)
    let p53 = db.proteins().find_by_gene("TP53")?;
    for p in &p53 {
        println!("Protein: {} | Length: {} aa", p.uniprot_accession, p.length);
        // ส่งต่อไปคำนวณด้วย scies-bio-th
        let bio_protein = scies_bio_th::Protein::new(p.sequence.as_bytes().to_vec());
    }

    // ค้นหาการกลายพันธุ์ทางคลินิก
    let variants = db.variants().find_by_gene("TP53")?;
    println!("Clinical variants in TP53: {}", variants.len());

    Ok(())
}
```

---

## 7. Dependency กับ `scies-bio-th`

| Type / Module จาก `scies-bio-th` | บทบาทใน `scies-bio-db` |
|---|---|
| `Dna` | Return type ของ `genome_query::fetch_sequence()` |
| `Protein` | สร้างจาก `ProteinEntry.sequence` เพื่อส่งให้ computation layer |
| `genome::Interval` | ใช้เป็น Input ของ Interval index |
| `variant::Variant` | ใช้เป็น base model สำหรับ `VariantRecord` |
| `annotation::Feature` | ใช้เป็น base model สำหรับ Gene / Transcript |
| `io::fasta` | Ingest pipeline อ่าน Reference genome |
| `io::gff3` | Ingest pipeline อ่าน Gene annotations |
| `io::vcf` | Ingest pipeline อ่าน Variant data |
| `proteomics::*` | คำนวณ Mass, Tryptic digest บน protein ที่ดึงมาจาก DB |

---

## 8. Roadmap

| Version | Scope | แหล่งข้อมูลที่รองรับ |
|---|---|---|
| **v0.1** | Taxonomy, Reference Genome (2-bit storage), Gene Models (GENCODE), Query API skeleton | NCBI RefSeq, Ensembl/GENCODE |
| **v0.2** | Proteomics (UniProt Human), Domain & GO Annotations | UniProtKB Swiss-Prot |
| **v0.3** | Variants (dbSNP, ClinVar), VariantEffect model | ClinVar, dbSNP |
| **v0.4** | Pathways (KEGG/Reactome), PPI (STRING), Gene Regulatory Network | KEGG, Reactome, STRING |
| **v0.5** | GWAS Catalog, Pharmacogenomics | GWAS Catalog, PharmGKB |
| **v0.6+** | Model Organisms: Mouse (10090), E. coli (562), Zebrafish (7955), Yeast (4932) | NCBI, Ensembl |

---

## 9. Out of Scope (v0.1)

- Genome Assembly pipeline
- Short-read mapping / alignment (อยู่ใน `scies-bio-th`)
- 3D Protein Structure (PDB)
- Real-time NCBI/Ensembl API calls (เก็บ local ก่อน)
- GPU/SIMD acceleration
