# The `.seq` 2-Bit Binary Sequence Format Specification (Version 2)

This document defines the formal binary format specification for `.seq` sequence files used by `scies-bio-db`.

---

## 1. Overview & Design Goals

The `.seq` binary format stores genomic sequences with maximal efficiency and strict integrity:
- **75% space reduction** over raw ASCII FASTA via 2-bit canonical nucleotide encoding.
- **Zero data loss**: Full 15-symbol IUPAC alphabet support (`A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N`) via a sorted ambiguity sidecar table.
- **Hardware-aligned 64-byte header** with magic signature `SCIESSEQ` and format version tagging.
- **Data integrity protection**: Built-in IEEE 802.3 CRC32 checksum over sequence and ambiguity payloads.
- **Crash-resilient atomic writes**: Multi-phase atomic write (`.tmp` + `sync_all` + atomic `rename`).
- **Constant-time $O(1)$ random-access slicing** with zero-copy `mmap` backing and $O(\log N_a)$ binary-search ambiguity recovery.
- **Backward compatibility**: Transparently decodes legacy unversioned v1 `.seq` files.

---

## 2. Binary File Layout (`.seq v2`)

A version 2 `.seq` file consists of a **64-byte fixed-size header**, followed by the **packed 2-bit canonical payload**, followed by the **sorted ambiguity sidecar payload**:

```text
+-------------------------------------------------------------------------+
| Offset (Bytes) | Field Name            | Type   | Description           |
+----------------+-----------------------+--------+-----------------------+
| 0..8           | Magic Bytes           | [u8;8] | ASCII "SCIESSEQ"      |
| 8..10          | Format Version        | u16 LE | 2                     |
| 10..12         | Flags                 | u16 LE | bit 0: has_ambiguity  |
|                |                       |        | bit 1: is_rna         |
| 12..20         | Sequence Length (N)   | u64 LE | Total base count      |
| 20..28         | Packed Length (L_p)   | u64 LE | ceil(N / 4) bytes     |
| 28..36         | Ambiguity Count (N_a) | u64 LE | Number of ambig bases |
| 36..44         | Ambiguity Length(L_a) | u64 LE | N_a * 9 bytes         |
| 44..48         | Payload CRC32         | u32 LE | IEEE 802.3 checksum   |
| 48..64         | Reserved / Padding    | [u8;16]| Zero-filled           |
+----------------+-----------------------+--------+-----------------------+
| 64..(64 + L_p) | Canonical 2-Bit Payload (ceil(N / 4) bytes)            |
|                | Packed MSB-first: 4 bases per byte                     |
+----------------+--------------------------------------------------------+
| (64 + L_p)..End| Ambiguity Sidecar Table (N_a * 9 bytes, sorted by pos) |
|                | Each record: [position: u64 LE (8B)] [iupac_code: u8]  |
+-------------------------------------------------------------------------+
```

### Total File Size Formula
$$\text{File Size} = 64 + \left\lceil \frac{N}{4} \right\rceil + (N_a \times 9) \text{ bytes}$$

Where:
- $N$ is total sequence length in bases ($0 \le N < 2^{64}$).
- $\lceil N / 4 \rceil$ is calculated as `(N + 3) / 4`.
- $N_a$ is the number of non-canonical IUPAC bases ($0 \le N_a \le N$).

---

## 3. Nucleotide Bit Encoding & IUPAC Sidecar

### Canonical 2-Bit Encoding
In the 2-bit packed sequence payload, each base is mapped to its closest canonical equivalent:

| Canonical Base | Bit Value | Integer | Notes |
|:---:|:---:|:---:|:---|
| **`A`** | `00` | `0` | Adenine |
| **`C`** | `01` | `1` | Cytosine |
| **`G`** | `10` | `2` | Guanine |
| **`T`** | `11` | `3` | Thymine |

### Lossless IUPAC Ambiguity Mapping
When non-canonical IUPAC characters or Uracil are present:
1. In the canonical 2-bit payload, the nearest match is stored so that downstream algorithms unaware of ambiguity can still operate reasonably:
   - `U` $\rightarrow$ `T` (`11`)
   - `R` (`A` or `G`) $\rightarrow$ `A` (`00`)
   - `Y` (`C` or `T`) $\rightarrow$ `C` (`01`)
   - `S` (`G` or `C`) $\rightarrow$ `G` (`10`)
   - `W` (`A` or `T`) $\rightarrow$ `T` (`11`)
   - `K` (`G` or `T`) $\rightarrow$ `G` (`10`)
   - `M` (`A` or `C`) $\rightarrow$ `A` (`00`)
   - `B` (`C`/`G`/`T`) $\rightarrow$ `C` (`01`)
   - `D` (`A`/`G`/`T`) $\rightarrow$ `A` (`00`)
   - `H` (`A`/`C`/`T`) $\rightarrow$ `A` (`00`)
   - `V` (`A`/`C`/`G`) $\rightarrow$ `A` (`00`)
   - `N` (any base) $\rightarrow$ `A` (`00`)
2. In the **ambiguity sidecar table**, the exact original ASCII character (e.g. `b'N'`, `b'U'`, `b'R'`) is recorded alongside its 0-based coordinate `(position: u64, code: u8)`.
3. Records in the ambiguity sidecar table are strictly sorted by `position`.

---

## 4. Slice Decoding & Ambiguity Recovery

When querying interval $[start, end)$:
1. Decode base characters from the 2-bit packed stream into an output buffer of size $end - start$.
2. If `N_a > 0`, locate the first ambiguity record where $position \ge start$ via **binary search** ($O(\log N_a)$).
3. Walk sequentially until $position \ge end$.
4. For each matching ambiguity record, overwrite `out[position - start] = code`.

This architecture ensures standard 4-base sequences incur zero overhead, while degenerate genomes retain **100% lossless fidelity**.

---

## 5. Integrity & Safety Guarantees

### IEEE 802.3 CRC32 Verification
- Checksum covers `packed_payload_bytes || ambiguity_payload_bytes`.
- Header verification ensures `mmap.len() >= 64 + L_p + L_a`.
- Checksum mismatch aborts with `BioDbError::ChecksumMismatch { expected, found }` before processing corrupt data.

### Atomic Writes
Writes follow a crash-resilient pattern:
1. Write to temporary file `<parent>/<chr>.seq.tmp.<pid>_<counter>`.
2. `writer.flush()` buffers to OS cache.
3. `file.sync_all()` persists data and metadata to durable physical storage.
4. `std::fs::rename` atomically swaps the temporary file into place.

### Path Traversal Protection
Path components (`assembly`, `chr`) are sanitized:
- Rejects paths containing `..`, `/`, `\`, or null bytes.
- Returns explicit `BioDbError::PathTraversal` error.

---

## 6. Legacy v1 Backward Compatibility

Existing files created under format version 1 (which begin with an 8-byte length header without magic bytes) are automatically detected and decoded using the legacy v1 decoder. All new writes automatically produce version 2 files.

