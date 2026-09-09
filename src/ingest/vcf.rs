//! VCF variant ingestion pipeline.
//!
//! Parses VCF records and inserts [`VariantRecord`] rows into the
//! [`SqliteStore`].

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use scies_bio_th::io::vcf::VcfReader;

use crate::model::variant::{VariantRecord, VariantType};
use crate::storage::SqliteStore;
use crate::Result;

/// Ingests genomic variants from VCF files into the [`SqliteStore`].
pub struct VcfIngestor;

impl VcfIngestor {
    /// Creates a new [`VcfIngestor`].
    pub fn new() -> Self {
        Self
    }

    /// Parses `path` as a VCF file and inserts all variant records into `store`
    /// scoped under `tax_id` and `assembly`.
    ///
    /// Returns the total number of variants successfully inserted.
    pub fn ingest_file(
        &self,
        path: &Path,
        tax_id: u32,
        assembly: &str,
        store: &SqliteStore,
    ) -> Result<u64> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let vcf_reader = VcfReader::new(reader);

        let mut count = 0u64;

        for item in vcf_reader {
            let record = item?;
            let ref_str = String::from_utf8_lossy(&record.reference).to_string();
            let alt_str = String::from_utf8_lossy(&record.alternate).to_string();
            let ref_len = record.reference.len();
            let alt_len = record.alternate.len();

            let vtype = if ref_len == 1 && alt_len == 1 {
                VariantType::Snp
            } else if ref_len < alt_len {
                VariantType::Insertion
            } else if ref_len > alt_len {
                VariantType::Deletion
            } else {
                VariantType::Mnp
            };

            let variant_id = record
                .id
                .filter(|s| s != "." && !s.is_empty())
                .unwrap_or_else(|| {
                    format!("{}:{}:{}>{}", record.chrom, record.pos, ref_str, alt_str)
                });

            let var_record = VariantRecord {
                variant_id,
                assembly_name: assembly.to_string(),
                chr: record.chrom.clone(),
                position: record.pos,
                reference: ref_str,
                alternate: alt_str,
                variant_type: vtype,
                tax_id,
            };

            store.insert_variant_record(&var_record)?;
            count += 1;
        }

        Ok(count)
    }
}

impl Default for VcfIngestor {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_vcf_ingestion() {
        let dir = tempdir().unwrap();
        let db_path = dir.path().join("test_vcf.db");
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

        let mut vcf_file = NamedTempFile::new().unwrap();
        writeln!(vcf_file, "##fileformat=VCFv4.2").unwrap();
        writeln!(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(vcf_file, "chr17\t7676154\trs1042522\tC\tG\t100\tPASS\t.").unwrap();
        vcf_file.flush().unwrap();

        let ingestor = VcfIngestor::new();
        let count = ingestor
            .ingest_file(vcf_file.path(), 9606, "GRCh38", &store)
            .unwrap();
        assert_eq!(count, 1);

        let var = store.find_variant_by_id("rs1042522").unwrap().unwrap();
        assert_eq!(var.chr, "chr17");
        // 1-based VCF POS 7676154 is converted to 0-based coordinate 7676153
        assert_eq!(var.position, 7676153);
        assert_eq!(var.reference, "C");
        assert_eq!(var.alternate, "G");
    }
}
