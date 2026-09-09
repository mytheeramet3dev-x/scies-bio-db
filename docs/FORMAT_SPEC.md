# The `.seq` 2-Bit Binary Sequence Format Specification

This document defines the formal binary format specification for `.seq` files used by `scies-bio-db`.

---

## 1. Overview

The `.seq` format is a compact, little-endian binary sequence serialization designed for:
- 75% space reduction compared to raw ASCII FASTA.
- Constant time $O(1)$ random access slicing without decompression or seeking overhead.
- Direct zero-copy memory mapping via the operating system page cache.

---

## 2. File Layout

A `.seq` file consists of an 8-byte fixed header followed by packed 2-bit nucleotide bytes:

```text
+-------------------------------------------------------------+
| Header (8 Bytes)                                            |
| Length of sequence in bases (N): 64-bit Little-Endian UInt  |
+-------------------------------------------------------------+
| Byte 8                                                      |
| Base[0] (Bits 7..6) | Base[1] (Bits 5..4)                   |
| Base[2] (Bits 3..2) | Base[3] (Bits 1..0)                   |
+-------------------------------------------------------------+
| Byte 9                                                      |
| Base[4] (Bits 7..6) | Base[5] (Bits 5..4)                   |
| Base[6] (Bits 3..2) | Base[7] (Bits 1..0)                   |
+-------------------------------------------------------------+
| ...                                                         |
+-------------------------------------------------------------+
| Byte 8 + ceil(N / 4) - 1 (Final Payload Byte)               |
| Contains remaining bases (with zero-padding in unused bits) |
+-------------------------------------------------------------+
```

### Total File Size Formula
$$\text{File Size} = 8 + \left\lceil \frac{N}{4} \right\rceil \text{ bytes}$$

Where:
- $N$ is the total sequence length in base pairs ($0 \le N < 2^{64}$).
- $\lceil N / 4 \rceil$ is calculated in integer arithmetic as `(N + 3) / 4`.

---

## 3. Nucleotide Bit Encoding Table

Each nucleotide is encoded into a 2-bit representation:

| Nucleotide Character | Bit Value (Binary) | Integer Value | Meaning |
|:---:|:---:|:---:|:---|
| **`A`** / **`a`** | `00` | `0` | Adenine |
| **`C`** / **`c`** | `01` | `1` | Cytosine |
| **`G`** / **`g`** | `10` | `2` | Guanine |
| **`T`** / **`t`** / **`U`** / **`u`** | `11` | `3` | Thymine / Uracil |

*Note*: Degenerate IUPAC characters (e.g. `N`, `R`, `Y`) or unknown bases are mapped deterministically to `A` (`00`) in the standard 2-bit mode.

---

## 4. Packing and Unpacking Algorithms

### Bit Packing (Encoding)
Bases are packed MSB-first (Most Significant Bit first) into bytes:

```text
Byte Layout:
Bit Index:   7   6   5   4   3   2   1   0
           [Base 0] [Base 1] [Base 2] [Base 3]
```

#### Encoding Logic in Rust:
```rust
pub fn encode(seq: &[u8]) -> Vec<u8> {
    let num_bytes = (seq.len() + 3) / 4;
    let mut out = Vec::with_capacity(num_bytes);

    for chunk in seq.chunks(4) {
        let mut byte = 0u8;
        for (i, &base) in chunk.iter().enumerate() {
            let bits = match base {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' | b'U' | b'u' => 0b11,
                _ => 0b00,
            };
            let shift = 6 - i * 2;
            byte |= bits << shift;
        }
        out.push(byte);
    }
    out
}
```

### Random Access Slice Decoding
To decode a half-open coordinate interval $[start, end)$ from memory-mapped data:

For any coordinate position $pos \in [start, end)$:
1. **Byte Index Calculation**:
   $$\text{byte\_idx} = \left\lfloor \frac{pos}{4} \right\rfloor$$
2. **Bit Shift Calculation**:
   $$\text{shift} = 6 - ((pos \bmod 4) \times 2)$$
3. **Bit Extraction**:
   $$\text{bits} = (\text{payload}[\text{byte\_idx}] \gg \text{shift}) \;\&\; 0\text{b}11$$

#### Decoding Logic in Rust:
```rust
pub fn decode(payload: &[u8], start: u64, end: u64) -> Vec<u8> {
    const BASES: [u8; 4] = *b"ACGT";
    let mut out = Vec::with_capacity((end - start) as usize);

    for pos in start..end {
        let byte_idx = (pos / 4) as usize;
        let shift = 6 - (pos % 4) * 2;
        let bits = (payload[byte_idx] >> shift) & 0b11;
        out.push(BASES[bits as usize]);
    }
    out
}
```

---

## 5. Padding & Boundary Conditions

When the total length $N$ is not divisible by 4:
- The last byte contains $N \bmod 4$ valid 2-bit pairs in the most significant bit positions.
- The remaining trailing bit positions in the final byte are zero-padded.
- When querying slices, any read beyond $N$ is guarded by the 8-byte length header and rejected with `BioDbError::InvalidCoordinate`.
