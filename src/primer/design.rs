use super::filter::{is_valid_primer_position, PrimerFilter, PrimerParams};
use super::thermodynamics::{calculate_gc_content, calculate_tm};
use crate::genome::fasta::GenomeReader;
use crate::genome::transcript::Transcript;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use tracing::debug;

/// A designed primer pair
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerPair {
    /// Gene symbol
    pub gene: String,
    /// Transcript accession
    pub transcript: String,
    /// Chromosome
    pub chromosome: String,
    /// Forward primer sequence (5' to 3')
    pub forward_seq: String,
    /// Reverse primer sequence (5' to 3')
    pub reverse_seq: String,
    /// Forward primer melting temperature
    pub forward_tm: f64,
    /// Reverse primer melting temperature
    pub reverse_tm: f64,
    /// Forward primer GC content (%)
    pub forward_gc: f64,
    /// Reverse primer GC content (%)
    pub reverse_gc: f64,
    /// Product size (amplicon length)
    pub product_size: usize,
    /// Number of introns spanned by the amplicon
    pub introns_spanned: usize,
    /// Genomic position of forward primer (5' end)
    pub forward_pos: u64,
    /// Genomic position of reverse primer (5' end)
    pub reverse_pos: u64,
    /// Exon index where forward primer starts (0-based in transcript coordinates)
    pub forward_exon: usize,
    /// Exon index where reverse primer starts (0-based in transcript coordinates)
    pub reverse_exon: usize,
}

/// Design primers for a transcript
///
/// Strategy:
/// 1. Get the spliced transcript sequence (concatenated exons)
/// 2. Find positions where primers span exon-exon junctions (for specificity)
/// 3. For each valid forward primer position, find valid reverse positions
/// 4. Filter primer pairs by quality criteria
pub fn design_primers(
    transcript: &Transcript,
    genome: &GenomeReader,
    params: &PrimerParams,
) -> Result<Vec<PrimerPair>> {
    let primer_len = params.primer_length;
    let filter = PrimerFilter::new(params.clone());

    // Get transcript sequence with position mapping
    let (seq, positions) = genome.get_transcript_sequence_with_positions(transcript)?;

    if seq.len() < primer_len * 2 + 50 {
        debug!(
            "Transcript {} too short ({} bp) for primer design",
            transcript.accession,
            seq.len()
        );
        return Ok(Vec::new());
    }

    // Build exon boundary map
    // Maps position in spliced sequence -> (exon_index, position_within_exon)
    let exon_map = build_exon_map(transcript, &positions);

    let mut primers = Vec::new();

    // Minimum product size (to span at least one intron effectively)
    let min_product = 80;
    let max_product = 300;

    // Find forward primer candidates
    for fwd_start in 0..seq.len().saturating_sub(primer_len + min_product) {
        let fwd_end = fwd_start + primer_len;
        let fwd_seq = &seq[fwd_start..fwd_end];

        // Check if this position is valid
        if !is_valid_primer_position(&seq, fwd_start, primer_len) {
            continue;
        }

        // Check if forward primer passes filter
        if !filter.passes(fwd_seq) {
            continue;
        }

        let fwd_tm = calculate_tm(fwd_seq);
        let fwd_gc = calculate_gc_content(fwd_seq);

        // Get exon info for forward primer
        let fwd_exon = exon_map.get(&fwd_start).map(|&(e, _)| e).unwrap_or(0);

        // Find reverse primer candidates
        let rev_search_start = fwd_end + min_product - primer_len;
        let rev_search_end = (fwd_end + max_product).min(seq.len());

        for rev_start in rev_search_start..rev_search_end.saturating_sub(primer_len) {
            let rev_end = rev_start + primer_len;

            if rev_end > seq.len() {
                break;
            }

            let rev_seq_fwd = &seq[rev_start..rev_end];

            // Check if position is valid
            if !is_valid_primer_position(&seq, rev_start, primer_len) {
                continue;
            }

            // Reverse primer needs to be reverse complemented
            let rev_seq = reverse_complement(rev_seq_fwd);

            // Check if reverse primer passes filter
            if !filter.passes(&rev_seq) {
                continue;
            }

            // Check primer pair compatibility (3' complementarity)
            if !filter.check_3prime_complementarity(fwd_seq, &rev_seq) {
                continue;
            }

            let rev_tm = calculate_tm(&rev_seq);
            let rev_gc = calculate_gc_content(&rev_seq);

            // Check Tm difference between primers
            if (fwd_tm - rev_tm).abs() > 5.0 {
                continue;
            }

            // Get exon info for reverse primer
            let rev_exon = exon_map.get(&rev_start).map(|&(e, _)| e).unwrap_or(0);

            // Calculate introns spanned
            let introns_spanned = if rev_exon > fwd_exon {
                rev_exon - fwd_exon
            } else {
                0
            };

            // We want primers that span at least one intron
            // Skip pairs where both primers are in the same exon
            if introns_spanned == 0 && transcript.exons.len() > 1 {
                continue;
            }

            let product_size = rev_end - fwd_start;

            // Get genomic positions
            let forward_pos = positions.get(fwd_start).copied().unwrap_or(0);
            let reverse_pos = positions.get(rev_start).copied().unwrap_or(0);

            primers.push(PrimerPair {
                gene: transcript.gene_symbol.clone(),
                transcript: transcript.accession.clone(),
                chromosome: transcript.chromosome.clone(),
                forward_seq: String::from_utf8_lossy(fwd_seq).to_string(),
                reverse_seq: String::from_utf8_lossy(&rev_seq).to_string(),
                forward_tm: round_2dp(fwd_tm),
                reverse_tm: round_2dp(rev_tm),
                forward_gc: round_2dp(fwd_gc),
                reverse_gc: round_2dp(rev_gc),
                product_size,
                introns_spanned,
                forward_pos,
                reverse_pos,
                forward_exon: fwd_exon,
                reverse_exon: rev_exon,
            });
        }
    }

    // Sort by product size (prefer smaller amplicons)
    primers.sort_by_key(|p| p.product_size);

    // Limit to best primers per transcript (avoid too many similar primers)
    let max_primers_per_transcript = 10;
    primers.truncate(max_primers_per_transcript);

    debug!(
        "Designed {} primers for transcript {}",
        primers.len(),
        transcript.accession
    );

    Ok(primers)
}

/// Build a map from spliced sequence position to (exon_index, position_within_exon)
fn build_exon_map(
    transcript: &Transcript,
    positions: &[u64],
) -> std::collections::HashMap<usize, (usize, usize)> {
    use std::collections::HashMap;

    let mut exon_map: HashMap<usize, (usize, usize)> = HashMap::new();

    // Sort exons by genomic position
    let mut sorted_exons: Vec<(usize, &crate::genome::transcript::Exon)> =
        transcript.exons.iter().enumerate().collect();
    sorted_exons.sort_by_key(|(_, e)| e.start);

    // For minus strand, reverse the exon order
    if transcript.strand == crate::genome::transcript::Strand::Minus {
        sorted_exons.reverse();
    }

    let mut seq_pos = 0;
    for (exon_idx, exon) in sorted_exons {
        let exon_len = (exon.end - exon.start) as usize;
        for i in 0..exon_len {
            exon_map.insert(seq_pos + i, (exon_idx, i));
        }
        seq_pos += exon_len;
    }

    exon_map
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

/// Round to 2 decimal places
fn round_2dp(val: f64) -> f64 {
    (val * 100.0).round() / 100.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT".to_vec());
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT".to_vec());
    }

    #[test]
    fn test_round_2dp() {
        assert_eq!(round_2dp(58.4287), 58.43);
        assert_eq!(round_2dp(58.421), 58.42);
    }
}
