# Tutorial & Production Simulation Pipelines

`scies-bio-db` ships with 4 complete, executable simulation pipelines located in the `examples/` directory.

---

## 🦁 Demo 1: Cave Lion De-Extinction & Artificial Womb
**File**: [`examples/cave_lion_deextinction_demo.rs`](../examples/cave_lion_deextinction_demo.rs)  
**Run Command**:
```bash
cargo run --example cave_lion_deextinction_demo
```

### Scientific Background
The Eurasian Cave Lion (*Panthera spelaea*, NCBI TaxID: `2562234`) was a dominant Pleistocene apex predator that went extinct approximately 14,000 years ago. Well-preserved permafrost mummies discovered in Yakutia, Siberia provide ancient somatic cell nuclei.

### The Pipeline Steps
```text
[Ancient Somatic Nucleus]
            │
            ▼
[IVG Synthetic Oocyte] ──► [Zona Pellucida Shell (ZP1-4)]
            │
            ▼
[SCNT Electrofusion] ──► [KDM4A Histone Demethylation (H3K9me3 barrier removed)]
            │
            ▼
[Pre-Implantation Culture] ──► [ZP Hatching at Day 6.0]
            │
            ▼
[Artificial Amniotic Biobag] ──► [110-Day Gestation] ──► [1.85 kg Live Neonate Cub]
```

1. **Taxonomy & Ancient Genome Registration**: Registers NCBI TaxID `2562234` (*Panthera spelaea*) into `meta.db`.
2. **Synthetic Oocyte Engineering (IVG)**: Generates a mature metaphase II oocyte with 100,000 maternal mitochondria and a synthetic **Zona Pellucida** protective glycoprotein matrix (thickness: 16.5 µm, elasticity: 18.5 kPa, containing ZP1, ZP2, ZP3, ZP4).
3. **SCNT & Epigenetic Barrier Erasure**: Simulates somatic nuclear transfer coupled with **KDM4A** mRNA microinjection, removing aberrant H3K9me3 marks and enabling Embryonic Genome Activation (EGA).
4. **Zona Pellucida Hatching**: Advances the embryo through cleavage (TwoCell $\rightarrow$ Morula $\rightarrow$ Blastocyst) to Day 6.0 hatching.
5. **Artificial Womb Extra-Uterine Gestation**: Translocates the hatched blastocyst into an artificial amniotic womb chamber (oxygenator $\text{PaO}_2 = 28\text{ mmHg}$, temperature $38.5^\circ\text{C}$). Gestates through gastrulation (Day 14), organogenesis (Day 30), fetal growth (Day 60), to a full-term 110-day healthy neonate cub weighing **1.85 kg**.
6. **Database Persistence**: Stores the de-extinction record into `seq/2562234/PanSpe_Yakutia_v1.0/ExUtero_Clone_Cub_01.seq`.

---

## 🎯 Demo 2: Personalized Cancer mRNA Neoantigen Vaccine Pipeline
**File**: [`examples/cancer_vaccine_demo.rs`](../examples/cancer_vaccine_demo.rs)  
**Run Command**:
```bash
cargo run --example cancer_vaccine_demo
```

### Scientific Background
Personalized mRNA vaccines instruct the patient's dendritic cells to present patient-specific tumor neoepitopes on MHC Class I molecules, priming cytotoxic CD8+ T-cells to eradicate cancer cells with minimal off-target autoimmunity.

### The Pipeline Steps
1. **Somatic Mutation Screening**: Scans key oncogenic driver mutations:
   - *BRAF* V600E (Valine $\rightarrow$ Glutamic Acid at position 600)
   - *TP53* R175H (Arginine $\rightarrow$ Histidine at position 175)
   - *KRAS* G12D (Glycine $\rightarrow$ Aspartate at position 12)
   - *PIK3CA* H1047R (Histidine $\rightarrow$ Arginine at position 1047)
2. **MHC Class I Binding Affinity Prediction**: Screens overlapping 9-mer peptides against patient HLA alleles (e.g. `HLA-A*02:01`). Computes IC50 (nM) and Differential Agretopic Index (DAI).
3. **Polytope Cassette Assembly**: Links top neoepitopes using proteasome-cleavable `AAY` linkers, flanked by a Secretion Signal peptide (SecSignal) and MHC Class I Trafficking Domain (MITD).
4. **Clinical mRNA Construct Synthesis**: Assembles clinical-grade mRNA:
   `[5' Cap][5' UTR + Kozak][Polytope Cassette][Dual Stop][Dual 3' UTR][Poly-A 120nt]`
   Optimizes human codons to achieve **69% GC content**.
5. **Catalog Archiving**: Saves the patient vaccine profile into `demo_biodb/seq/9606/` for downstream clinical tracking.

---

## 🥚 Demo 3: In Vitro Gametogenesis (IVG) & SCNT
**File**: [`examples/ivg_scnt_demo.rs`](../examples/ivg_scnt_demo.rs)  
**Run Command**:
```bash
cargo run --example ivg_scnt_demo
```

### Scientific Background
In Vitro Gametogenesis (IVG) differentiates pluripotent stem cells into primordial germ cell-like cells (PGCLCs) and functional gametes (oocytes or spermatozoa), enabling reproduction independent of adult donor availability.

---

## 🧬 Demo 4: Authentic Multi-Species DNA Ingestion
**File**: [`examples/real_dna_demo.rs`](../examples/real_dna_demo.rs)  
**Run Command**:
```bash
cargo run --example real_dna_demo
```

### Verified Sequences
Ingests authentic reference sequences from NCBI GenBank:
- **Human Insulin (*INS*)**: $333\text{ bp}$ coding region on chromosome 11.
- **Human Beta-Globin (*HBB*)**: $444\text{ bp}$ coding region on chromosome 11.
- **E. coli $\beta$-Galactosidase (*lacZ*)**: $3,075\text{ bp}$ on *E. coli* K-12 str. MG1655.
- **SARS-CoV-2 Spike Glycoprotein (RBD)**: $2,835\text{ bp}$ from Wuhan-Hu-1 (`NC_045512.2`).

Verifies exact byte-level retrieval, GC content calculation, and demonstrates disk compression efficiency (~3.8x to 4x smaller than raw ASCII text).
