use serde::{Deserialize, Serialize};

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            '+' => Some(Strand::Plus),
            '-' => Some(Strand::Minus),
            _ => None,
        }
    }

    pub fn as_char(&self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
        }
    }
}

/// An exon within a transcript
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Exon {
    /// 0-based start position (inclusive)
    pub start: u64,
    /// 0-based end position (exclusive)
    pub end: u64,
    /// Exon number (1-based, from 5' to 3' of transcript)
    pub number: usize,
}

impl Exon {
    pub fn new(start: u64, end: u64, number: usize) -> Self {
        Self { start, end, number }
    }

    /// Length of the exon in base pairs
    pub fn length(&self) -> u64 {
        self.end - self.start
    }
}

/// A transcript (single isoform of a gene)
#[derive(Debug, Clone)]
pub struct Transcript {
    /// Transcript accession (e.g., "NM_000546")
    pub accession: String,
    /// Gene symbol (e.g., "TP53")
    pub gene_symbol: String,
    /// Chromosome (e.g., "chr17")
    pub chromosome: String,
    /// Strand orientation
    pub strand: Strand,
    /// Transcription start (0-based, inclusive)
    pub tx_start: u64,
    /// Transcription end (0-based, exclusive)
    pub tx_end: u64,
    /// Coding sequence start (0-based, inclusive)
    pub cds_start: u64,
    /// Coding sequence end (0-based, exclusive)
    pub cds_end: u64,
    /// Exons ordered by genomic position
    pub exons: Vec<Exon>,
}

impl Transcript {
    /// Check if this transcript is coding (has a valid CDS)
    pub fn is_coding(&self) -> bool {
        self.cds_start < self.cds_end
    }

    /// Get exons in 5' to 3' order (accounts for strand)
    pub fn exons_5to3(&self) -> Vec<&Exon> {
        let mut exons: Vec<&Exon> = self.exons.iter().collect();
        match self.strand {
            Strand::Plus => exons.sort_by_key(|e| e.start),
            Strand::Minus => exons.sort_by_key(|e| std::cmp::Reverse(e.start)),
        }
        exons
    }

    /// Get total exonic length
    pub fn exonic_length(&self) -> u64 {
        self.exons.iter().map(|e| e.length()).sum()
    }

    /// Get the number of introns
    pub fn intron_count(&self) -> usize {
        if self.exons.is_empty() {
            0
        } else {
            self.exons.len() - 1
        }
    }

    /// Calculate how many introns are between two genomic positions
    /// Positions should be on the same chromosome
    pub fn introns_between(&self, pos1: u64, pos2: u64) -> usize {
        let (start, end) = if pos1 < pos2 {
            (pos1, pos2)
        } else {
            (pos2, pos1)
        };

        // Count exon boundaries between the two positions
        let mut intron_count = 0;
        let sorted_exons: Vec<&Exon> = {
            let mut e: Vec<&Exon> = self.exons.iter().collect();
            e.sort_by_key(|ex| ex.start);
            e
        };

        for i in 0..sorted_exons.len().saturating_sub(1) {
            let intron_start = sorted_exons[i].end;
            let intron_end = sorted_exons[i + 1].start;

            // Check if this intron is fully contained between our positions
            if intron_start > start && intron_end < end {
                intron_count += 1;
            }
        }

        intron_count
    }
}

/// A gene with all its transcripts
#[derive(Debug, Clone)]
pub struct Gene {
    /// Gene symbol (e.g., "TP53")
    pub symbol: String,
    /// All transcripts for this gene
    pub transcripts: Vec<Transcript>,
}

impl Gene {
    pub fn new(symbol: String) -> Self {
        Self {
            symbol,
            transcripts: Vec::new(),
        }
    }

    /// Add a transcript to this gene
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.push(transcript);
    }

    /// Get the canonical transcript (longest coding sequence)
    pub fn canonical_transcript(&self) -> Option<&Transcript> {
        self.transcripts
            .iter()
            .filter(|t| t.is_coding())
            .max_by_key(|t| t.cds_end - t.cds_start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exon_length() {
        let exon = Exon::new(100, 200, 1);
        assert_eq!(exon.length(), 100);
    }

    #[test]
    fn test_strand_conversion() {
        assert_eq!(Strand::from_char('+'), Some(Strand::Plus));
        assert_eq!(Strand::from_char('-'), Some(Strand::Minus));
        assert_eq!(Strand::from_char('x'), None);
        assert_eq!(Strand::Plus.as_char(), '+');
        assert_eq!(Strand::Minus.as_char(), '-');
    }

    #[test]
    fn test_transcript_introns() {
        let transcript = Transcript {
            accession: "NM_TEST".to_string(),
            gene_symbol: "TEST".to_string(),
            chromosome: "chr1".to_string(),
            strand: Strand::Plus,
            tx_start: 100,
            tx_end: 500,
            cds_start: 150,
            cds_end: 450,
            exons: vec![
                Exon::new(100, 150, 1),
                Exon::new(200, 250, 2),
                Exon::new(300, 350, 3),
                Exon::new(400, 500, 4),
            ],
        };

        assert_eq!(transcript.intron_count(), 3);
        // Between positions 125 and 425, we cross 2 full introns (150-200 and 250-300)
        // The intron 350-400 is also fully contained
        assert_eq!(transcript.introns_between(125, 425), 2);
    }
}
