//! primer-genie: A Rust CLI tool for designing gene-specific primers
//!
//! This library provides functionality for:
//! - Downloading reference genomes and annotations
//! - Parsing FASTA and refGene annotation files
//! - Designing PCR primers that span intron-exon boundaries
//! - Checking primer specificity across the transcriptome
//! - Outputting primers in JSONL format

pub mod commands;
pub mod config;
pub mod genome;
pub mod output;
pub mod primer;

// Re-export commonly used types
pub use config::{AssemblyConfig, Species};
pub use genome::{Gene, GenomeReader, Transcript};
pub use output::PrimerOutput;
pub use primer::{PrimerFilter, PrimerPair, PrimerParams};
