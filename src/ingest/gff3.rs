//! GFF3 annotation ingestion pipeline.
//!
//! Parses GFF3 and GTF records and inserts gene, transcript, and exon rows into the
//! [`SqliteStore`].

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use scies_bio_th::io::gff::GffReader;

use crate::model::gene::{Biotype, Exon, Gene, GeneFeature, Transcript, TranscriptBiotype};
use crate::model::genome::{Assembly, Strand};
use crate::storage::SqliteStore;
use crate::Result;

/// Detailed report of a GFF3 / GTF annotation ingestion (Section 5 DB3).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Gff3IngestReport {
    /// Number of gene features successfully ingested.
    pub genes_count: u64,
    /// Number of transcript features successfully ingested.
    pub transcripts_count: u64,
    /// Number of exon features successfully ingested.
    pub exons_count: u64,
    /// Number of malformed or incomplete features skipped.
    pub skipped_count: u64,
    /// Total features inserted across all types.
    pub total_inserted: u64,
}

/// Ingests genomic feature annotations from GFF3 / GTF files into the [`SqliteStore`].
pub struct Gff3Ingestor;

impl Gff3Ingestor {
    /// Creates a new [`Gff3Ingestor`].
    pub fn new() -> Self {
        Self
    }

    /// Parses `path` as a GFF3 / GTF file and inserts all recognised features into `store`
    /// scoped under `tax_id` and `assembly`.
    ///
    /// Converts GFF3 1-based inclusive interval [start, end] into 0-based half-open [start - 1, end).
    ///
    /// Returns the total number of features (genes, transcripts, exons) successfully inserted.
    pub fn ingest_file(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
        store: &SqliteStore,
    ) -> Result<u64> {
        let report = self.ingest_file_report(path, tax_id, assembly, store)?;
        Ok(report.total_inserted)
    }

    /// Parses `path` as a GFF3 / GTF file and returns a detailed [`Gff3IngestReport`]
    /// reporting counts of ingested genes, transcripts, exons, and skipped records.
    pub fn ingest_file_report(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
        store: &SqliteStore,
    ) -> Result<Gff3IngestReport> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let gff_reader = GffReader::new(reader);

        // Auto-register assembly in SQLite assemblies catalog if not already present
        let asm = Assembly {
            tax_id,
            assembly_name: assembly.to_string(),
            source: Some(path.to_string_lossy().to_string()),
            created_at: None,
            notes: None,
        };
        let _ = store.insert_assembly(&asm);

        let mut report = Gff3IngestReport::default();

        let mut gene_features = Vec::new();
        let mut genes = Vec::new();
        let mut transcripts = Vec::new();
        let mut exons = Vec::new();
        let mut transcript_exon_counts: HashMap<String, u32> = HashMap::new();

        for item in gff_reader {
            let record = match item {
                Ok(r) => r,
                Err(_) => {
                    report.skipped_count += 1;
                    continue;
                }
            };

            // Coordinate integrity check: start must be <= end
            if record.start > record.end {
                report.skipped_count += 1;
                continue;
            }

            let ftype = record.feature_type.to_ascii_lowercase();
            let attrs = parse_attributes(&record.attributes);

            let strand = match record.strand {
                scies_bio_th::genome::Strand::Reverse => Strand::Reverse,
                _ => Strand::Forward,
            };

            let source_str = if record.source.is_empty() {
                None
            } else {
                Some(record.source.clone())
            };

            match ftype.as_str() {
                "gene" | "pseudogene" => {
                    let gene_id = attrs
                        .get("id")
                        .or_else(|| attrs.get("gene_id"))
                        .cloned()
                        .unwrap_or_else(|| {
                            format!("{}_{}_{}", record.seqid, record.start, record.end)
                        });
                    let clean_id = gene_id
                        .strip_prefix("gene:")
                        .unwrap_or(&gene_id)
                        .to_string();

                    // Section 5 DB3: support symbol from gene_name, Name, gene, symbol
                    let symbol = attrs
                        .get("gene_name")
                        .or_else(|| attrs.get("name"))
                        .or_else(|| attrs.get("gene"))
                        .or_else(|| attrs.get("symbol"))
                        .cloned()
                        .unwrap_or_else(|| clean_id.clone());

                    let name = attrs
                        .get("description")
                        .or_else(|| attrs.get("name"))
                        .cloned()
                        .unwrap_or_else(|| symbol.clone());

                    let biotype_str = attrs
                        .get("biotype")
                        .or_else(|| attrs.get("gene_biotype"))
                        .map(|s| s.as_str())
                        .unwrap_or(if ftype == "pseudogene" {
                            "pseudogene"
                        } else {
                            "protein_coding"
                        });

                    let biotype = match biotype_str {
                        "protein_coding" => Biotype::ProteinCoding,
                        "lncrna" | "lnc_rna" => Biotype::LongNonCodingRna,
                        "mirna" => Biotype::MicroRna,
                        "snrna" => Biotype::SmallNuclearRna,
                        "pseudogene" => Biotype::Pseudogene,
                        "rrna" => Biotype::RibosomalRna,
                        other => Biotype::Other(other.to_string()),
                    };

                    gene_features.push(GeneFeature {
                        tax_id,
                        assembly_name: assembly.to_string(),
                        gene_id: clean_id.clone(),
                        symbol: symbol.clone(),
                        chrom: record.seqid.clone(),
                        start_bp: record.start,
                        end_bp: record.end,
                        strand,
                        feature_type: "gene".to_string(),
                        source: source_str,
                    });

                    genes.push(Gene {
                        gene_id: clean_id,
                        assembly_name: assembly.to_string(),
                        symbol,
                        name,
                        biotype,
                        chr: record.seqid,
                        start: record.start,
                        end: record.end,
                        strand,
                        tax_id,
                    });
                }
                "mrna" | "transcript" | "lnc_rna" | "mirna" | "ncrna" => {
                    let transcript_id = attrs
                        .get("id")
                        .or_else(|| attrs.get("transcript_id"))
                        .cloned()
                        .unwrap_or_else(|| {
                            format!("{}_{}_{}", record.seqid, record.start, record.end)
                        });
                    let clean_tid = transcript_id
                        .strip_prefix("transcript:")
                        .unwrap_or(&transcript_id)
                        .to_string();

                    let raw_parent = attrs
                        .get("parent")
                        .or_else(|| attrs.get("gene_id"))
                        .cloned()
                        .unwrap_or_default();
                    let gene_id = raw_parent
                        .strip_prefix("gene:")
                        .unwrap_or(&raw_parent)
                        .to_string();

                    let biotype = match ftype.as_str() {
                        "mrna" => TranscriptBiotype::MRNA,
                        "lnc_rna" | "mirna" | "ncrna" => TranscriptBiotype::NCRNA,
                        "pseudogene" => TranscriptBiotype::Pseudogene,
                        other => TranscriptBiotype::Other(other.to_string()),
                    };

                    transcripts.push(Transcript {
                        transcript_id: clean_tid,
                        gene_id,
                        assembly_name: assembly.to_string(),
                        biotype,
                        chr: record.seqid,
                        start: record.start,
                        end: record.end,
                        strand,
                        exon_count: 0,
                    });
                }
                "exon" => {
                    let exon_id = attrs
                        .get("id")
                        .or_else(|| attrs.get("exon_id"))
                        .cloned()
                        .unwrap_or_else(|| {
                            format!("exon_{}_{}_{}", record.seqid, record.start, record.end)
                        });
                    let clean_eid = exon_id
                        .strip_prefix("exon:")
                        .unwrap_or(&exon_id)
                        .to_string();

                    let raw_parent = attrs
                        .get("parent")
                        .or_else(|| attrs.get("transcript_id"))
                        .cloned()
                        .unwrap_or_default();
                    let transcript_id = raw_parent
                        .strip_prefix("transcript:")
                        .unwrap_or(&raw_parent)
                        .to_string();

                    let exon_number = attrs
                        .get("exon_number")
                        .or_else(|| attrs.get("rank"))
                        .and_then(|s| s.parse::<u32>().ok())
                        .unwrap_or_else(|| {
                            let count = transcript_exon_counts
                                .entry(transcript_id.clone())
                                .or_insert(0);
                            *count += 1;
                            *count
                        });

                    exons.push(Exon {
                        exon_id: clean_eid,
                        transcript_id,
                        assembly_name: assembly.to_string(),
                        chr: record.seqid,
                        start: record.start,
                        end: record.end,
                        strand,
                        exon_number,
                    });
                }
                _ => {
                    // Non-target features (e.g. chromosome, region, biological_region)
                }
            }
        }

        // Fast O(E + T) aggregation of exon counts per transcript
        let mut exon_counts: HashMap<&str, u32> = HashMap::with_capacity(transcripts.len());
        for exon in &exons {
            *exon_counts.entry(exon.transcript_id.as_str()).or_insert(0) += 1;
        }
        for t in &mut transcripts {
            t.exon_count = exon_counts
                .get(t.transcript_id.as_str())
                .copied()
                .unwrap_or(1);
        }

        // Batch insert GeneFeature into genes catalog
        store.insert_gene_features_batch(&gene_features)?;
        report.genes_count = gene_features.len() as u64;

        // Best-effort batch insert into transcripts and exons tables
        if let Err(_) = store.insert_transcripts_batch(&transcripts) {
            // If foreign keys fail for transcripts due to missing parent genes, insert individually
            for t in &transcripts {
                if store.insert_transcript(t).is_ok() {
                    report.transcripts_count += 1;
                } else {
                    report.skipped_count += 1;
                }
            }
        } else {
            report.transcripts_count = transcripts.len() as u64;
        }

        if let Err(_) = store.insert_exons_batch(&exons) {
            for e in &exons {
                if store.insert_exon(e).is_ok() {
                    report.exons_count += 1;
                } else {
                    report.skipped_count += 1;
                }
            }
        } else {
            report.exons_count = exons.len() as u64;
        }

        report.total_inserted = report.genes_count + report.transcripts_count + report.exons_count;
        Ok(report)
    }
}

impl Default for Gff3Ingestor {
    fn default() -> Self {
        Self::new()
    }
}

/// Helper to parse GFF3 and GTF attribute columns (e.g. `ID=gene1;Name=BRCA1` or `gene_id "ENSG1"; gene_name "BRCA1";`).
fn parse_attributes(attrs: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in attrs.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((k, v)) = part.split_once('=') {
            let val = v.trim().trim_matches('"');
            map.insert(k.trim().to_ascii_lowercase(), val.to_string());
        } else if let Some((k, v)) = part.split_once(' ') {
            let val = v.trim().trim_matches('"');
            map.insert(k.trim().to_ascii_lowercase(), val.to_string());
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_gff3_ingestion() {
        let dir = tempdir().unwrap();
        let db_path = dir.path().join("test_gff3.db");
        let store = SqliteStore::open(&db_path).unwrap();
        store.init_schema().unwrap();

        // Register organism and genome first
        let org = crate::model::taxonomy::Organism {
            tax_id: 9606,
            scientific_name: "Homo sapiens".to_string(),
            common_name: Some("Human".to_string()),
            lineage: crate::model::taxonomy::Lineage {
                domain: "Eukaryota".into(),
                kingdom: "Animalia".into(),
                phylum: "Chordata".into(),
                class: "Mammalia".into(),
                order: "Primates".into(),
                family: "Hominidae".into(),
                genus: "Homo".into(),
                species: "sapiens".into(),
            },
        };
        store.insert_organism(&org).unwrap();

        let genome = crate::model::genome::ReferenceGenome {
            assembly_name: "GRCh38".to_string(),
            tax_id: 9606,
            status: crate::model::genome::AssemblyStatus::Primary,
            release_year: 2013,
        };
        store.insert_reference_genome(&genome).unwrap();

        let mut gff_file = NamedTempFile::new().unwrap();
        writeln!(gff_file, "##gff-version 3").unwrap();
        writeln!(
            gff_file,
            "chr17\tEnsembl\tgene\t7668402\t7687550\t.\t+\t.\tID=gene:ENSG00000141510;Name=TP53;biotype=protein_coding"
        )
        .unwrap();
        writeln!(
            gff_file,
            "chr17\tEnsembl\tmRNA\t7668402\t7687550\t.\t+\t.\tID=transcript:ENST00000269305;Parent=gene:ENSG00000141510;biotype=protein_coding"
        )
        .unwrap();
        writeln!(
            gff_file,
            "chr17\tEnsembl\texon\t7668402\t7669600\t.\t+\t.\tID=exon:ENSE000001;Parent=transcript:ENST00000269305;exon_number=1"
        )
        .unwrap();
        gff_file.flush().unwrap();

        let ingestor = Gff3Ingestor::new();
        let count = ingestor
            .ingest_file(gff_file.path(), 9606, "GRCh38", &store)
            .unwrap();
        assert_eq!(count, 3);

        // Verify gene query
        let gene = store.find_gene_by_symbol(9606, "TP53").unwrap().unwrap();
        assert_eq!(gene.gene_id, "ENSG00000141510");
        assert_eq!(gene.symbol, "TP53");
        // Coordinate conversion: 1-based 7668402 -> 0-based 7668401
        assert_eq!(gene.start, 7668401);
        assert_eq!(gene.end, 7687550);

        // Verify transcript
        let transcripts = store.find_transcripts_by_gene("ENSG00000141510").unwrap();
        assert_eq!(transcripts.len(), 1);
        assert_eq!(transcripts[0].transcript_id, "ENST00000269305");
        assert_eq!(transcripts[0].exon_count, 1);
    }

    #[test]
    fn test_ingest_mini_genes_gff3_fixture() {
        let dir = tempdir().unwrap();
        let db_path = dir.path().join("test_mini_genes.db");
        let store = SqliteStore::open(&db_path).unwrap();
        store.init_schema().unwrap();

        let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let gff_path = manifest_dir.join("testdata/mini_genes.gff3");

        let ingestor = Gff3Ingestor::new();
        let rep = ingestor
            .ingest_file_report(&gff_path, 9606, "T2T-CHM13v2.0", &store)
            .unwrap();

        assert_eq!(rep.genes_count, 3);
        assert_eq!(rep.transcripts_count, 3);
        assert_eq!(rep.exons_count, 3);
        assert_eq!(rep.skipped_count, 1); // The bad coordinate line (5000 > 4000) was skipped
        assert_eq!(rep.total_inserted, 9);

        // 1. Check TERT (gene_name=TERT, forward strand)
        let tert = store
            .get_gene_feature_by_symbol(9606, "T2T-CHM13v2.0", "TERT")
            .unwrap()
            .expect("TERT must be found");
        assert_eq!(tert.symbol, "TERT");
        assert_eq!(tert.chrom, "chr5");
        assert_eq!(tert.start_bp, 1294800);
        assert_eq!(tert.end_bp, 1295600);
        assert_eq!(tert.strand, Strand::Forward);

        // 2. Check SIRT6 (Name=SIRT6, forward strand)
        let sirt6 = store
            .get_gene_feature_by_symbol(9606, "T2T-CHM13v2.0", "SIRT6")
            .unwrap()
            .expect("SIRT6 must be found");
        assert_eq!(sirt6.symbol, "SIRT6");
        assert_eq!(sirt6.chrom, "chr19");
        assert_eq!(sirt6.start_bp, 4182000);
        assert_eq!(sirt6.end_bp, 4182800);
        assert_eq!(sirt6.strand, Strand::Forward);

        // 3. Check TP53 (gene=TP53, reverse strand)
        let tp53 = store
            .get_gene_feature_by_symbol(9606, "T2T-CHM13v2.0", "TP53")
            .unwrap()
            .expect("TP53 must be found");
        assert_eq!(tp53.symbol, "TP53");
        assert_eq!(tp53.chrom, "chr17");
        assert_eq!(tp53.start_bp, 7687000);
        assert_eq!(tp53.end_bp, 7688000);
        assert_eq!(tp53.strand, Strand::Reverse);

        // 4. Region query test
        let in_region = store
            .list_gene_features_in_region(9606, "T2T-CHM13v2.0", "chr5", 1294000, 1296000)
            .unwrap();
        assert_eq!(in_region.len(), 1);
        assert_eq!(in_region[0].symbol, "TERT");
    }
}
