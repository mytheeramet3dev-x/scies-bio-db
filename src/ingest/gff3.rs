//! GFF3 annotation ingestion pipeline.
//!
//! Parses GFF3 records and inserts gene, transcript, and exon rows into the
//! [`SqliteStore`].

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use scies_bio_th::io::gff::GffReader;

use crate::model::gene::{Biotype, Exon, Gene, Transcript, TranscriptBiotype};
use crate::model::genome::Strand;
use crate::storage::SqliteStore;
use crate::Result;

/// Ingests genomic feature annotations from GFF3 files into the [`SqliteStore`].
pub struct Gff3Ingestor;

impl Gff3Ingestor {
    /// Creates a new [`Gff3Ingestor`].
    pub fn new() -> Self {
        Self
    }

    /// Parses `path` as a GFF3 file and inserts all recognised features into `store`
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
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let gff_reader = GffReader::new(reader);

        let mut inserted_count = 0u64;

        let mut genes = Vec::new();
        let mut transcripts = Vec::new();
        let mut exons = Vec::new();
        let mut transcript_exon_counts: HashMap<String, u32> = HashMap::new();

        for item in gff_reader {
            let record = item?;
            let ftype = record.feature_type.to_ascii_lowercase();
            let attrs = parse_attributes(&record.attributes);

            let strand = match record.strand {
                scies_bio_th::genome::Strand::Reverse => Strand::Reverse,
                _ => Strand::Forward,
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

                    let symbol = attrs
                        .get("name")
                        .or_else(|| attrs.get("gene_name"))
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
                _ => {}
            }
        }

        // Update transcript exon counts
        for t in &mut transcripts {
            let count = exons
                .iter()
                .filter(|e| e.transcript_id == t.transcript_id)
                .count() as u32;
            t.exon_count = if count > 0 { count } else { 1 };
        }

        // Insert into database
        for gene in &genes {
            store.insert_gene(gene)?;
            inserted_count += 1;
        }
        for transcript in &transcripts {
            store.insert_transcript(transcript)?;
            inserted_count += 1;
        }
        for exon in &exons {
            store.insert_exon(exon)?;
            inserted_count += 1;
        }

        Ok(inserted_count)
    }
}

impl Default for Gff3Ingestor {
    fn default() -> Self {
        Self::new()
    }
}

/// Helper to parse GFF3 attribute column (e.g. `ID=gene1;Name=BRCA1`).
fn parse_attributes(attrs: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in attrs.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        if let Some((k, v)) = part.split_once('=') {
            map.insert(k.trim().to_ascii_lowercase(), v.trim().to_string());
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
}
