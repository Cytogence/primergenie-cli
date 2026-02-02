use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use tracing::{debug, warn};

use super::transcript::{Exon, Gene, Strand, Transcript};

/// Parse a refGene.txt file and return a map of gene symbol -> Gene
///
/// refGene.txt format (tab-delimited):
/// bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames
///
/// - name = transcript accession (NM_*)
/// - name2 = gene symbol
pub fn parse_refgene(path: &Path) -> Result<HashMap<String, Gene>> {
    let file = File::open(path).context("Failed to open refGene file")?;

    let reader: Box<dyn BufRead> = if path
        .extension()
        .map_or(false, |ext| ext == "gz")
    {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut genes: HashMap<String, Gene> = HashMap::new();
    let mut line_num = 0;
    let mut skipped = 0;

    for line in reader.lines() {
        line_num += 1;
        let line = line.context("Failed to read line")?;

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        match parse_refgene_line(&line) {
            Ok(transcript) => {
                let gene_symbol = transcript.gene_symbol.clone();
                genes
                    .entry(gene_symbol.clone())
                    .or_insert_with(|| Gene::new(gene_symbol))
                    .add_transcript(transcript);
            }
            Err(e) => {
                debug!("Skipping line {}: {}", line_num, e);
                skipped += 1;
            }
        }
    }

    if skipped > 0 {
        warn!("Skipped {} malformed lines in refGene file", skipped);
    }

    debug!(
        "Parsed {} genes with {} transcripts",
        genes.len(),
        genes.values().map(|g| g.transcripts.len()).sum::<usize>()
    );

    Ok(genes)
}

/// Parse a single line from refGene.txt
fn parse_refgene_line(line: &str) -> Result<Transcript> {
    let fields: Vec<&str> = line.split('\t').collect();

    if fields.len() < 16 {
        anyhow::bail!("Expected at least 16 fields, got {}", fields.len());
    }

    // Field indices (0-based):
    // 0: bin, 1: name (accession), 2: chrom, 3: strand, 4: txStart, 5: txEnd
    // 6: cdsStart, 7: cdsEnd, 8: exonCount, 9: exonStarts, 10: exonEnds
    // 11: score, 12: name2 (gene symbol), 13: cdsStartStat, 14: cdsEndStat, 15: exonFrames

    let accession = fields[1].to_string();
    let chromosome = fields[2].to_string();
    let strand = Strand::from_char(fields[3].chars().next().unwrap_or('+'))
        .ok_or_else(|| anyhow::anyhow!("Invalid strand: {}", fields[3]))?;
    let tx_start: u64 = fields[4].parse().context("Invalid txStart")?;
    let tx_end: u64 = fields[5].parse().context("Invalid txEnd")?;
    let cds_start: u64 = fields[6].parse().context("Invalid cdsStart")?;
    let cds_end: u64 = fields[7].parse().context("Invalid cdsEnd")?;
    let exon_count: usize = fields[8].parse().context("Invalid exonCount")?;
    let gene_symbol = fields[12].to_string();

    // Parse exon starts and ends (comma-separated, with trailing comma)
    let exon_starts: Vec<u64> = fields[9]
        .trim_end_matches(',')
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<u64>())
        .collect::<Result<Vec<_>, _>>()
        .context("Invalid exonStarts")?;

    let exon_ends: Vec<u64> = fields[10]
        .trim_end_matches(',')
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<u64>())
        .collect::<Result<Vec<_>, _>>()
        .context("Invalid exonEnds")?;

    if exon_starts.len() != exon_count || exon_ends.len() != exon_count {
        anyhow::bail!(
            "Exon count mismatch: expected {}, got {} starts and {} ends",
            exon_count,
            exon_starts.len(),
            exon_ends.len()
        );
    }

    // Build exons
    let mut exons: Vec<Exon> = exon_starts
        .iter()
        .zip(exon_ends.iter())
        .enumerate()
        .map(|(i, (&start, &end))| Exon::new(start, end, i + 1))
        .collect();

    // Sort by position
    exons.sort_by_key(|e| e.start);

    // Renumber based on strand (5' to 3')
    if strand == Strand::Minus {
        let n = exons.len();
        for (i, exon) in exons.iter_mut().enumerate() {
            exon.number = n - i;
        }
    }

    Ok(Transcript {
        accession,
        gene_symbol,
        chromosome,
        strand,
        tx_start,
        tx_end,
        cds_start,
        cds_end,
        exons,
    })
}

/// Filter genes to only include those in the provided list
pub fn filter_genes(genes: HashMap<String, Gene>, gene_list: &[String]) -> HashMap<String, Gene> {
    let gene_set: std::collections::HashSet<&str> =
        gene_list.iter().map(|s| s.as_str()).collect();

    genes
        .into_iter()
        .filter(|(symbol, _)| gene_set.contains(symbol.as_str()))
        .collect()
}

/// Read a gene list file (one gene symbol per line)
pub fn read_gene_list(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path).context("Failed to open gene list file")?;
    let reader = BufReader::new(file);

    let genes: Vec<String> = reader
        .lines()
        .filter_map(|line| {
            line.ok().and_then(|l| {
                let trimmed = l.trim().to_string();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    None
                } else {
                    Some(trimmed)
                }
            })
        })
        .collect();

    Ok(genes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_refgene_line() {
        // Example refGene line
        let line = "585\tNM_001005484\tchr1\t+\t69090\t70008\t69090\t70008\t1\t69090,\t70008,\t0\tOR4F5\tcmpl\tcmpl\t0,";
        let transcript = parse_refgene_line(line).unwrap();

        assert_eq!(transcript.accession, "NM_001005484");
        assert_eq!(transcript.gene_symbol, "OR4F5");
        assert_eq!(transcript.chromosome, "chr1");
        assert_eq!(transcript.strand, Strand::Plus);
        assert_eq!(transcript.tx_start, 69090);
        assert_eq!(transcript.tx_end, 70008);
        assert_eq!(transcript.exons.len(), 1);
    }

    #[test]
    fn test_parse_refgene_line_multiple_exons() {
        // Multi-exon gene
        let line = "0\tNM_000546\tchr17\t-\t7668402\t7687550\t7669609\t7676594\t11\t7668402,7670609,7673534,7673700,7674180,7674858,7675052,7675993,7676381,7676520,7687376,\t7669690,7670715,7673608,7673837,7674290,7674971,7675236,7676272,7676403,7676622,7687550,\t0\tTP53\tcmpl\tcmpl\t2,2,1,0,0,0,2,1,1,1,-1,";
        let transcript = parse_refgene_line(line).unwrap();

        assert_eq!(transcript.accession, "NM_000546");
        assert_eq!(transcript.gene_symbol, "TP53");
        assert_eq!(transcript.strand, Strand::Minus);
        assert_eq!(transcript.exons.len(), 11);
        assert_eq!(transcript.intron_count(), 10);
    }
}
