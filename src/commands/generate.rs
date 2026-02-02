use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;
use tracing::{info, warn};

use crate::config::{self, Species};
use crate::genome::annotations::{filter_genes, parse_refgene, read_gene_list};
use crate::genome::fasta::GenomeReader;
use crate::genome::transcript::{Gene, Transcript};
use crate::output::jsonl::JsonlWriter;
use crate::output::PrimerOutput;
use crate::primer::{design_primers, PrimerPair, PrimerParams, SpecificityChecker};

/// Parameters for the generate command
pub struct GenerateParams {
    pub species: String,
    pub assembly: Option<String>,
    pub output: Option<String>,
    pub primer_length: usize,
    pub min_tm: f64,
    pub max_tm: f64,
    pub min_gc: f64,
    pub max_gc: f64,
    pub max_homopolymer: usize,
    pub genes_file: Option<String>,
    pub threads: Option<usize>,
    pub skip_specificity: bool,
}

/// Run the generate command
pub async fn run(params: GenerateParams) -> Result<()> {
    let species = Species::from_str(&params.species)?;
    let assembly = params
        .assembly
        .as_deref()
        .unwrap_or(species.default_assembly());

    info!(
        "Generating primers for {} ({})",
        species.name(),
        assembly
    );

    // Check that reference files exist
    if !config::references_exist(species, assembly)? {
        anyhow::bail!(
            "Reference files not found. Run 'primer-genie fetch {}' first.",
            params.species
        );
    }

    // Set up thread pool
    if let Some(threads) = params.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }

    // Build primer parameters
    let primer_params = PrimerParams {
        primer_length: params.primer_length,
        min_tm: params.min_tm,
        max_tm: params.max_tm,
        min_gc: params.min_gc,
        max_gc: params.max_gc,
        max_homopolymer: params.max_homopolymer,
        ..PrimerParams::default()
    };

    // Load annotations
    info!("Loading gene annotations...");
    let annotations_path = config::annotations_path(species, assembly)?;
    let mut genes = parse_refgene(&annotations_path)?;
    info!("Loaded {} genes", genes.len());

    // Filter genes if a gene list was provided
    if let Some(genes_file) = &params.genes_file {
        let gene_list = read_gene_list(&PathBuf::from(genes_file))?;
        info!("Filtering to {} genes from list", gene_list.len());
        genes = filter_genes(genes, &gene_list);
        info!("After filtering: {} genes", genes.len());
    }

    if genes.is_empty() {
        warn!("No genes to process");
        return Ok(());
    }

    // Load genome
    info!("Loading genome (this may take a while for large genomes)...");
    let genome_path = config::genome_path(species, assembly)?;
    let genome = GenomeReader::load(&genome_path)?;

    // Phase 1: Design primers for all transcripts
    info!("Phase 1: Designing primers...");
    let all_primers = design_all_primers(&genes, &genome, &primer_params)?;
    info!("Designed {} primer pairs", all_primers.len());

    if all_primers.is_empty() {
        warn!("No primers designed");
        return Ok(());
    }

    // Phase 2: Specificity checking
    let primers_with_targets = if params.skip_specificity {
        info!("Skipping specificity check (--skip-specificity)");
        // Set all targets to 1 (unknown specificity)
        all_primers
            .into_iter()
            .map(|p| (p, 1usize))
            .collect::<Vec<_>>()
    } else {
        info!("Phase 2: Checking primer specificity...");
        check_specificity(&genes, &genome, all_primers)?
    };

    // Phase 3: Output
    info!("Phase 3: Writing output...");
    let output_path = params.output.as_ref().map(|s| PathBuf::from(s));
    write_output(
        primers_with_targets,
        species.name(),
        assembly,
        output_path.as_deref(),
    )?;

    info!("Done!");

    Ok(())
}

/// Design primers for all genes in parallel
fn design_all_primers(
    genes: &HashMap<String, Gene>,
    genome: &GenomeReader,
    params: &PrimerParams,
) -> Result<Vec<PrimerPair>> {
    let genes_vec: Vec<&Gene> = genes.values().collect();

    let pb = ProgressBar::new(genes_vec.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} genes ({eta})")?
            .progress_chars("#>-"),
    );

    let all_primers = Mutex::new(Vec::new());

    genes_vec.par_iter().for_each(|gene| {
        // Process each transcript for this gene
        for transcript in &gene.transcripts {
            // Skip non-coding transcripts
            if !transcript.is_coding() {
                continue;
            }

            // Skip transcripts on unrecognized chromosomes
            if !genome.has_chromosome(&transcript.chromosome) {
                continue;
            }

            match design_primers(transcript, genome, params) {
                Ok(primers) => {
                    if !primers.is_empty() {
                        let mut all = all_primers.lock().unwrap();
                        all.extend(primers);
                    }
                }
                Err(e) => {
                    tracing::debug!(
                        "Failed to design primers for {}: {}",
                        transcript.accession,
                        e
                    );
                }
            }
        }

        pb.inc(1);
    });

    pb.finish_with_message("Primer design complete");

    let primers = all_primers.into_inner().unwrap();
    Ok(primers)
}

/// Check specificity for all designed primers
fn check_specificity(
    genes: &HashMap<String, Gene>,
    genome: &GenomeReader,
    primers: Vec<PrimerPair>,
) -> Result<Vec<(PrimerPair, usize)>> {
    // Collect all transcripts for the specificity index
    let all_transcripts: Vec<&Transcript> = genes
        .values()
        .flat_map(|g| g.transcripts.iter())
        .filter(|t| t.is_coding() && genome.has_chromosome(&t.chromosome))
        .collect();

    // Build specificity checker
    let checker = SpecificityChecker::new(&all_transcripts, genome)?;

    // Check specificity in batches
    let pb = ProgressBar::new(primers.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} primers ({eta})")?
            .progress_chars("#>-"),
    );

    let batch_size = 1000;
    let mut results: Vec<(PrimerPair, usize)> = Vec::with_capacity(primers.len());

    for chunk in primers.chunks(batch_size) {
        let targets = checker.check_batch(chunk, 500); // Max 500bp product for specificity

        for (primer, target_count) in chunk.iter().zip(targets.iter()) {
            results.push((primer.clone(), *target_count));
        }

        pb.inc(chunk.len() as u64);
    }

    pb.finish_with_message("Specificity check complete");

    // Log statistics
    let specific = results.iter().filter(|(_, t)| *t == 1).count();
    let multi_target = results.iter().filter(|(_, t)| *t > 1).count();
    info!(
        "Specificity results: {} specific (1 target), {} multi-target (>1 target)",
        specific, multi_target
    );

    Ok(results)
}

/// Write primers to output
fn write_output(
    primers: Vec<(PrimerPair, usize)>,
    species: &str,
    assembly: &str,
    output_path: Option<&std::path::Path>,
) -> Result<()> {
    let mut writer = JsonlWriter::new(output_path)?;

    for (primer, targets) in primers {
        let output = PrimerOutput::from_primer_pair(&primer, species, assembly, targets);
        writer.write(&output)?;
    }

    writer.flush()?;

    let path_display = output_path
        .map(|p| p.display().to_string())
        .unwrap_or_else(|| "stdout".to_string());

    info!("Output written to {}", path_display);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primer_params_default() {
        let params = PrimerParams::default();
        assert_eq!(params.primer_length, 20);
        assert_eq!(params.min_tm, 50.0);
        assert_eq!(params.max_tm, 65.0);
    }
}
