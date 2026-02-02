use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use tracing::{debug, info};

use super::transcript::{Strand, Transcript};

/// Index entry for a chromosome/contig in the FASTA file
#[derive(Debug, Clone)]
pub struct FastaIndexEntry {
    /// Sequence name (e.g., "chr1")
    pub name: String,
    /// Total length in bases
    pub length: u64,
    /// Byte offset to first base in file
    pub offset: u64,
    /// Number of bases per line
    pub line_bases: u64,
    /// Number of bytes per line (including newline)
    pub line_width: u64,
}

/// Reader for indexed FASTA files
/// Supports both gzipped and uncompressed FASTA files
pub struct GenomeReader {
    /// Chromosome sequences loaded into memory
    sequences: HashMap<String, Vec<u8>>,
}

impl GenomeReader {
    /// Load a genome FASTA file into memory
    ///
    /// For gzipped files, the entire file is decompressed into memory.
    /// This is necessary because gzip doesn't support random access.
    pub fn load(path: &Path) -> Result<Self> {
        info!("Loading genome from {}...", path.display());

        let file = File::open(path).context("Failed to open genome file")?;

        let reader: Box<dyn BufRead> = if path
            .extension()
            .map_or(false, |ext| ext == "gz")
        {
            info!("Decompressing gzipped genome (this may take a while)...");
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
        let mut current_name: Option<String> = None;
        let mut current_seq: Vec<u8> = Vec::new();

        for line in reader.lines() {
            let line = line.context("Failed to read line")?;

            if line.starts_with('>') {
                // Save previous sequence
                if let Some(name) = current_name.take() {
                    debug!("Loaded {} ({} bp)", name, current_seq.len());
                    sequences.insert(name, std::mem::take(&mut current_seq));
                }

                // Start new sequence
                let name = line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                current_name = Some(name);
                current_seq = Vec::with_capacity(250_000_000); // Pre-allocate for large chromosomes
            } else {
                // Sequence data - convert to uppercase bytes
                current_seq.extend(line.trim().as_bytes().iter().map(|b| b.to_ascii_uppercase()));
            }
        }

        // Save last sequence
        if let Some(name) = current_name {
            debug!("Loaded {} ({} bp)", name, current_seq.len());
            sequences.insert(name, current_seq);
        }

        info!("Loaded {} chromosomes/contigs", sequences.len());

        Ok(Self { sequences })
    }

    /// Get a subsequence from the genome
    ///
    /// - `chrom`: Chromosome name (e.g., "chr1")
    /// - `start`: 0-based start position (inclusive)
    /// - `end`: 0-based end position (exclusive)
    ///
    /// Returns the sequence as uppercase ASCII bytes
    pub fn get_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let seq = self
            .sequences
            .get(chrom)
            .ok_or_else(|| anyhow::anyhow!("Chromosome not found: {}", chrom))?;

        let start = start as usize;
        let end = end as usize;

        if start >= seq.len() || end > seq.len() || start >= end {
            anyhow::bail!(
                "Invalid coordinates for {}: {}..{} (length: {})",
                chrom,
                start,
                end,
                seq.len()
            );
        }

        Ok(seq[start..end].to_vec())
    }

    /// Get the reverse complement of a sequence
    pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'G' | b'g' => b'C',
                b'C' | b'c' => b'G',
                b'N' | b'n' => b'N',
                _ => b'N',
            })
            .collect()
    }

    /// Extract the spliced transcript sequence (concatenated exons)
    ///
    /// Returns the mRNA sequence in 5' to 3' orientation
    pub fn get_transcript_sequence(&self, transcript: &Transcript) -> Result<Vec<u8>> {
        let mut sequence = Vec::new();

        // Get exons in 5' to 3' order
        let exons = transcript.exons_5to3();

        for exon in exons {
            let exon_seq = self.get_sequence(&transcript.chromosome, exon.start, exon.end)?;
            sequence.extend(exon_seq);
        }

        // Reverse complement if on minus strand
        if transcript.strand == Strand::Minus {
            sequence = Self::reverse_complement(&sequence);
        }

        Ok(sequence)
    }

    /// Get the concatenated exon sequence with genomic position mapping
    ///
    /// Returns:
    /// - The spliced sequence (5' to 3')
    /// - A vector mapping each position in the spliced sequence to its genomic coordinate
    pub fn get_transcript_sequence_with_positions(
        &self,
        transcript: &Transcript,
    ) -> Result<(Vec<u8>, Vec<u64>)> {
        let mut sequence = Vec::new();
        let mut positions = Vec::new();

        // Get exons in genomic order first
        let mut sorted_exons: Vec<_> = transcript.exons.iter().collect();
        sorted_exons.sort_by_key(|e| e.start);

        for exon in &sorted_exons {
            let exon_seq = self.get_sequence(&transcript.chromosome, exon.start, exon.end)?;
            sequence.extend(exon_seq);

            // Map each base to its genomic position
            for pos in exon.start..exon.end {
                positions.push(pos);
            }
        }

        // Reverse if on minus strand
        if transcript.strand == Strand::Minus {
            sequence = Self::reverse_complement(&sequence);
            positions.reverse();
        }

        Ok((sequence, positions))
    }

    /// Check if a chromosome exists in the genome
    pub fn has_chromosome(&self, chrom: &str) -> bool {
        self.sequences.contains_key(chrom)
    }

    /// Get the length of a chromosome
    pub fn chromosome_length(&self, chrom: &str) -> Option<u64> {
        self.sequences.get(chrom).map(|s| s.len() as u64)
    }

    /// Get list of all chromosome names
    pub fn chromosomes(&self) -> Vec<&str> {
        self.sequences.keys().map(|s| s.as_str()).collect()
    }
}

/// Parse a FASTA index (.fai) file
pub fn parse_fasta_index(path: &Path) -> Result<HashMap<String, FastaIndexEntry>> {
    let file = File::open(path).context("Failed to open FASTA index")?;
    let reader = BufReader::new(file);

    let mut index = HashMap::new();

    for line in reader.lines() {
        let line = line.context("Failed to read line")?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 5 {
            continue;
        }

        let entry = FastaIndexEntry {
            name: fields[0].to_string(),
            length: fields[1].parse().context("Invalid length")?,
            offset: fields[2].parse().context("Invalid offset")?,
            line_bases: fields[3].parse().context("Invalid line_bases")?,
            line_width: fields[4].parse().context("Invalid line_width")?,
        };

        index.insert(entry.name.clone(), entry);
    }

    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(
            GenomeReader::reverse_complement(b"ATCG"),
            b"CGAT".to_vec()
        );
        assert_eq!(
            GenomeReader::reverse_complement(b"AAAA"),
            b"TTTT".to_vec()
        );
        assert_eq!(
            GenomeReader::reverse_complement(b"GCGC"),
            b"GCGC".to_vec()
        );
    }
}
