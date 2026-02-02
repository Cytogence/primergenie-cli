use anyhow::Result;
use clap::{Parser, Subcommand};
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

mod commands;
mod config;
mod genome;
mod output;
mod primer;

use commands::{fetch, generate};

#[derive(Parser)]
#[command(name = "primer-genie")]
#[command(author, version, about = "A Rust CLI tool for designing gene-specific primers")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Download reference files for a species
    Fetch {
        /// Species to fetch (human or mouse)
        species: String,

        /// Genome assembly version (e.g., hg38, hg19, mm39, mm10)
        #[arg(long, short = 'a')]
        assembly: Option<String>,
    },

    /// Generate primers for a species
    Generate {
        /// Species to generate primers for (human or mouse)
        species: String,

        /// Genome assembly version (e.g., hg38, hg19, mm39, mm10)
        #[arg(long, short = 'a')]
        assembly: Option<String>,

        /// Output file (default: stdout)
        #[arg(long, short = 'o')]
        output: Option<String>,

        /// Primer length in base pairs
        #[arg(long, default_value = "20")]
        primer_length: usize,

        /// Minimum melting temperature (°C)
        #[arg(long, default_value = "50.0")]
        min_tm: f64,

        /// Maximum melting temperature (°C)
        #[arg(long, default_value = "65.0")]
        max_tm: f64,

        /// Minimum GC content (%)
        #[arg(long, default_value = "50.0")]
        min_gc: f64,

        /// Maximum GC content (%)
        #[arg(long, default_value = "65.0")]
        max_gc: f64,

        /// Maximum consecutive identical bases (homopolymer)
        #[arg(long, default_value = "3")]
        max_homopolymer: usize,

        /// Only process genes listed in this file (one per line)
        #[arg(long)]
        genes: Option<String>,

        /// Number of worker threads
        #[arg(long, short = 't')]
        threads: Option<usize>,

        /// Skip specificity check (faster, but includes multi-target primers)
        #[arg(long)]
        skip_specificity: bool,
    },
}

#[tokio::main]
async fn main() -> Result<()> {
    // Initialize logging
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(EnvFilter::from_default_env().add_directive(tracing::Level::INFO.into()))
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Fetch { species, assembly } => {
            fetch::run(&species, assembly.as_deref()).await?;
        }
        Commands::Generate {
            species,
            assembly,
            output,
            primer_length,
            min_tm,
            max_tm,
            min_gc,
            max_gc,
            max_homopolymer,
            genes,
            threads,
            skip_specificity,
        } => {
            let params = generate::GenerateParams {
                species,
                assembly,
                output,
                primer_length,
                min_tm,
                max_tm,
                min_gc,
                max_gc,
                max_homopolymer,
                genes_file: genes,
                threads,
                skip_specificity,
            };
            generate::run(params).await?;
        }
    }

    Ok(())
}
