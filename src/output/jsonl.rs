use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::primer::PrimerPair;

/// Output format for primers in JSONL
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerOutput {
    /// Gene symbol
    pub gene: String,
    /// Transcript accession
    pub transcript: String,
    /// Species name
    pub species: String,
    /// Genome assembly
    pub assembly: String,
    /// Chromosome
    pub chromosome: String,
    /// Forward primer sequence (5' to 3')
    pub forward_seq: String,
    /// Reverse primer sequence (5' to 3')
    pub reverse_seq: String,
    /// Forward primer melting temperature (°C)
    pub forward_tm: f64,
    /// Reverse primer melting temperature (°C)
    pub reverse_tm: f64,
    /// Forward primer GC content (%)
    pub forward_gc: f64,
    /// Reverse primer GC content (%)
    pub reverse_gc: f64,
    /// Amplicon size (bp)
    pub product_size: usize,
    /// Number of introns spanned by the amplicon
    pub introns_spanned: usize,
    /// Genomic position of forward primer (5' end)
    pub forward_pos: u64,
    /// Genomic position of reverse primer (5' end)
    pub reverse_pos: u64,
    /// Number of genomic locations this primer pair matches (1 = specific)
    pub targets: usize,
}

impl PrimerOutput {
    /// Create a PrimerOutput from a PrimerPair
    pub fn from_primer_pair(
        pair: &PrimerPair,
        species: &str,
        assembly: &str,
        targets: usize,
    ) -> Self {
        Self {
            gene: pair.gene.clone(),
            transcript: pair.transcript.clone(),
            species: species.to_string(),
            assembly: assembly.to_string(),
            chromosome: pair.chromosome.clone(),
            forward_seq: pair.forward_seq.clone(),
            reverse_seq: pair.reverse_seq.clone(),
            forward_tm: pair.forward_tm,
            reverse_tm: pair.reverse_tm,
            forward_gc: pair.forward_gc,
            reverse_gc: pair.reverse_gc,
            product_size: pair.product_size,
            introns_spanned: pair.introns_spanned,
            forward_pos: pair.forward_pos,
            reverse_pos: pair.reverse_pos,
            targets,
        }
    }
}

/// Write primers to a JSONL file
pub fn write_primers_jsonl(primers: &[PrimerOutput], path: &Path) -> Result<()> {
    let file = File::create(path).context("Failed to create output file")?;
    let mut writer = BufWriter::new(file);

    for primer in primers {
        let json = serde_json::to_string(primer).context("Failed to serialize primer")?;
        writeln!(writer, "{}", json).context("Failed to write primer")?;
    }

    writer.flush().context("Failed to flush output")?;

    Ok(())
}

/// Write primers to stdout in JSONL format
pub fn write_primers_jsonl_stdout(primers: &[PrimerOutput]) -> Result<()> {
    let stdout = std::io::stdout();
    let mut writer = BufWriter::new(stdout.lock());

    for primer in primers {
        let json = serde_json::to_string(primer).context("Failed to serialize primer")?;
        writeln!(writer, "{}", json).context("Failed to write primer")?;
    }

    writer.flush().context("Failed to flush output")?;

    Ok(())
}

/// Streaming writer for JSONL output (for memory efficiency with large datasets)
pub struct JsonlWriter {
    writer: BufWriter<Box<dyn Write + Send>>,
}

impl JsonlWriter {
    /// Create a new JSONL writer
    pub fn new(path: Option<&Path>) -> Result<Self> {
        let writer: Box<dyn Write + Send> = match path {
            Some(p) => {
                let file = File::create(p).context("Failed to create output file")?;
                Box::new(file)
            }
            None => Box::new(std::io::stdout()),
        };

        Ok(Self {
            writer: BufWriter::new(writer),
        })
    }

    /// Write a single primer to the output
    pub fn write(&mut self, primer: &PrimerOutput) -> Result<()> {
        let json = serde_json::to_string(primer).context("Failed to serialize primer")?;
        writeln!(self.writer, "{}", json).context("Failed to write primer")?;
        Ok(())
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("Failed to flush output")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primer::PrimerPair;

    #[test]
    fn test_primer_output_serialization() {
        let pair = PrimerPair {
            gene: "TP53".to_string(),
            transcript: "NM_000546".to_string(),
            chromosome: "chr17".to_string(),
            forward_seq: "GCTATCGGATCCGCAGCCCC".to_string(),
            reverse_seq: "TGAAGAACATCGCTGTAGAT".to_string(),
            forward_tm: 58.42,
            reverse_tm: 59.18,
            forward_gc: 65.0,
            reverse_gc: 55.0,
            product_size: 85,
            introns_spanned: 1,
            forward_pos: 7675000,
            reverse_pos: 7675085,
            forward_exon: 0,
            reverse_exon: 1,
        };

        let output = PrimerOutput::from_primer_pair(&pair, "human", "hg38", 1);

        let json = serde_json::to_string(&output).unwrap();
        assert!(json.contains("\"gene\":\"TP53\""));
        assert!(json.contains("\"species\":\"human\""));
        assert!(json.contains("\"assembly\":\"hg38\""));
        assert!(json.contains("\"targets\":1"));
    }
}
