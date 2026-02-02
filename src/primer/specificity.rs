use aho_corasick::AhoCorasick;
use anyhow::Result;
use std::collections::HashMap;
use tracing::{debug, info};

use super::PrimerPair;
use crate::genome::fasta::GenomeReader;
use crate::genome::transcript::Transcript;

/// Specificity checker for primer pairs
///
/// Uses Aho-Corasick algorithm for efficient multi-pattern matching
/// to count how many transcripts each primer pair could amplify.
pub struct SpecificityChecker {
    /// Map from transcript accession to its sequence
    transcript_sequences: HashMap<String, Vec<u8>>,
}

impl SpecificityChecker {
    /// Build a specificity checker from a set of transcripts
    pub fn new(
        transcripts: &[&Transcript],
        genome: &GenomeReader,
    ) -> Result<Self> {
        info!("Building transcript sequence index for specificity checking...");

        let mut transcript_sequences = HashMap::new();

        for transcript in transcripts {
            match genome.get_transcript_sequence(transcript) {
                Ok(seq) => {
                    transcript_sequences.insert(transcript.accession.clone(), seq);
                }
                Err(e) => {
                    debug!(
                        "Skipping transcript {} for specificity index: {}",
                        transcript.accession, e
                    );
                }
            }
        }

        info!(
            "Built index with {} transcript sequences",
            transcript_sequences.len()
        );

        Ok(Self {
            transcript_sequences,
        })
    }

    /// Count how many transcripts a primer pair could amplify
    ///
    /// A primer pair matches a transcript if:
    /// 1. The forward primer matches somewhere in the transcript (5' to 3')
    /// 2. The reverse primer matches downstream of the forward (as reverse complement)
    /// 3. The distance between them is reasonable (< max_product_size)
    pub fn count_targets(&self, primer_pair: &PrimerPair, max_product_size: usize) -> usize {
        let forward_seq = primer_pair.forward_seq.as_bytes();
        let reverse_seq_rc = reverse_complement(primer_pair.reverse_seq.as_bytes());

        let mut target_count = 0;

        for (_accession, transcript_seq) in &self.transcript_sequences {
            if self.primer_matches_transcript(
                forward_seq,
                &reverse_seq_rc,
                transcript_seq,
                max_product_size,
            ) {
                target_count += 1;
            }
        }

        target_count
    }

    /// Check if a primer pair could amplify a given transcript
    fn primer_matches_transcript(
        &self,
        forward: &[u8],
        reverse_rc: &[u8],
        transcript: &[u8],
        max_product_size: usize,
    ) -> bool {
        // Find all forward primer matches
        let forward_matches = find_matches(forward, transcript);

        if forward_matches.is_empty() {
            return false;
        }

        // Find all reverse primer matches (already reverse complemented)
        let reverse_matches = find_matches(reverse_rc, transcript);

        if reverse_matches.is_empty() {
            return false;
        }

        // Check if any forward/reverse pair has valid spacing
        for &fwd_pos in &forward_matches {
            for &rev_pos in &reverse_matches {
                // Reverse primer should be downstream of forward
                if rev_pos > fwd_pos {
                    let distance = rev_pos - fwd_pos + reverse_rc.len();
                    if distance <= max_product_size {
                        return true;
                    }
                }
            }
        }

        false
    }

    /// Check specificity for multiple primer pairs efficiently using Aho-Corasick
    pub fn check_batch(
        &self,
        primer_pairs: &[PrimerPair],
        max_product_size: usize,
    ) -> Vec<usize> {
        if primer_pairs.is_empty() {
            return Vec::new();
        }

        // Build patterns for Aho-Corasick
        // Each primer pair contributes 2 patterns: forward and reverse_rc
        let mut patterns: Vec<Vec<u8>> = Vec::with_capacity(primer_pairs.len() * 2);
        let mut pattern_to_pair: Vec<(usize, bool)> = Vec::with_capacity(primer_pairs.len() * 2);

        for (pair_idx, pair) in primer_pairs.iter().enumerate() {
            // Forward primer
            patterns.push(pair.forward_seq.as_bytes().to_vec());
            pattern_to_pair.push((pair_idx, true)); // true = forward

            // Reverse primer (reverse complement)
            patterns.push(reverse_complement(pair.reverse_seq.as_bytes()));
            pattern_to_pair.push((pair_idx, false)); // false = reverse
        }

        let ac = match AhoCorasick::new(&patterns) {
            Ok(ac) => ac,
            Err(_) => {
                // Fall back to individual checking
                return primer_pairs
                    .iter()
                    .map(|p| self.count_targets(p, max_product_size))
                    .collect();
            }
        };

        // For each primer pair, track (transcript, fwd_positions, rev_positions)
        let mut pair_matches: Vec<HashMap<String, (Vec<usize>, Vec<usize>)>> =
            vec![HashMap::new(); primer_pairs.len()];

        // Scan all transcripts
        for (accession, seq) in &self.transcript_sequences {
            // Find all pattern matches in this transcript
            for mat in ac.find_iter(seq) {
                let pattern_idx = mat.pattern().as_usize();
                let (pair_idx, is_forward) = pattern_to_pair[pattern_idx];
                let pos = mat.start();

                let entry = pair_matches[pair_idx]
                    .entry(accession.clone())
                    .or_insert_with(|| (Vec::new(), Vec::new()));

                if is_forward {
                    entry.0.push(pos);
                } else {
                    entry.1.push(pos);
                }
            }
        }

        // Count valid targets for each primer pair
        let mut target_counts = vec![0usize; primer_pairs.len()];

        for (pair_idx, matches) in pair_matches.iter().enumerate() {
            let rev_len = primer_pairs[pair_idx].reverse_seq.len();

            for (_accession, (fwd_positions, rev_positions)) in matches {
                if fwd_positions.is_empty() || rev_positions.is_empty() {
                    continue;
                }

                // Check if any combination has valid spacing
                let mut is_valid = false;
                'outer: for &fwd_pos in fwd_positions {
                    for &rev_pos in rev_positions {
                        if rev_pos > fwd_pos {
                            let distance = rev_pos - fwd_pos + rev_len;
                            if distance <= max_product_size {
                                is_valid = true;
                                break 'outer;
                            }
                        }
                    }
                }

                if is_valid {
                    target_counts[pair_idx] += 1;
                }
            }
        }

        target_counts
    }
}

/// Find all positions where pattern matches in text
fn find_matches(pattern: &[u8], text: &[u8]) -> Vec<usize> {
    let mut positions = Vec::new();

    if pattern.is_empty() || text.len() < pattern.len() {
        return positions;
    }

    for i in 0..=text.len() - pattern.len() {
        if &text[i..i + pattern.len()] == pattern {
            positions.push(i);
        }
    }

    positions
}

/// Reverse complement a DNA sequence
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_matches() {
        let text = b"ATCGATCGATCG";
        let pattern = b"ATCG";

        let matches = find_matches(pattern, text);
        assert_eq!(matches, vec![0, 4, 8]);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT".to_vec());
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT".to_vec());
    }
}
