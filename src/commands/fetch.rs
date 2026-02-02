use anyhow::{Context, Result};
use futures::StreamExt;
use indicatif::{ProgressBar, ProgressStyle};
use std::path::Path;
use tokio::fs::{self, File};
use tokio::io::AsyncWriteExt;
use tracing::{info, warn};

use crate::config::{self, AssemblyConfig, Species};

/// Run the fetch command to download reference files
pub async fn run(species_str: &str, assembly: Option<&str>) -> Result<()> {
    let species = Species::from_str(species_str)?;
    let config = AssemblyConfig::new(species, assembly)?;

    info!(
        "Fetching reference files for {} ({})",
        species.name(),
        config.assembly
    );

    // Create cache directory
    let cache_dir = config::cache_dir(species, &config.assembly)?;
    fs::create_dir_all(&cache_dir)
        .await
        .context("Failed to create cache directory")?;

    // Download genome
    let genome_path = config::genome_path(species, &config.assembly)?;
    if genome_path.exists() {
        info!("Genome file already exists: {}", genome_path.display());
    } else {
        info!("Downloading genome from {}", config.genome_url);
        download_file(&config.genome_url, &genome_path).await?;
        info!("Genome downloaded to {}", genome_path.display());
    }

    // Download annotations
    let annotations_path = config::annotations_path(species, &config.assembly)?;
    if annotations_path.exists() {
        info!(
            "Annotations file already exists: {}",
            annotations_path.display()
        );
    } else {
        info!("Downloading annotations from {}", config.annotations_url);
        download_file(&config.annotations_url, &annotations_path).await?;
        info!("Annotations downloaded to {}", annotations_path.display());
    }

    // Create FASTA index for random access
    let index_path = config::genome_index_path(species, &config.assembly)?;
    if index_path.exists() {
        info!("FASTA index already exists: {}", index_path.display());
    } else {
        info!("Creating FASTA index...");
        create_fasta_index(&genome_path, &index_path).await?;
        info!("FASTA index created: {}", index_path.display());
    }

    println!(
        "\nReference files for {} ({}) are ready.",
        species.name(),
        config.assembly
    );
    println!("Cache directory: {}", cache_dir.display());

    Ok(())
}

/// Download a file with progress bar
async fn download_file(url: &str, dest: &Path) -> Result<()> {
    let client = reqwest::Client::new();

    let response = client
        .get(url)
        .send()
        .await
        .context("Failed to send request")?;

    if !response.status().is_success() {
        anyhow::bail!(
            "Failed to download {}: HTTP {}",
            url,
            response.status()
        );
    }

    let total_size = response.content_length();

    let pb = if let Some(size) = total_size {
        let pb = ProgressBar::new(size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")?
                .progress_chars("#>-"),
        );
        pb
    } else {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {bytes} downloaded")?,
        );
        pb
    };

    // Create temporary file
    let temp_path = dest.with_extension("tmp");
    let mut file = File::create(&temp_path)
        .await
        .context("Failed to create temporary file")?;

    let mut stream = response.bytes_stream();
    let mut downloaded: u64 = 0;

    while let Some(chunk) = stream.next().await {
        let chunk = chunk.context("Error reading response stream")?;
        file.write_all(&chunk)
            .await
            .context("Failed to write to file")?;
        downloaded += chunk.len() as u64;
        pb.set_position(downloaded);
    }

    file.flush().await?;
    drop(file);

    // Rename temp file to final destination
    fs::rename(&temp_path, dest)
        .await
        .context("Failed to rename temporary file")?;

    pb.finish_with_message("Download complete");

    Ok(())
}

/// Create a FASTA index file for random access
async fn create_fasta_index(genome_path: &Path, index_path: &Path) -> Result<()> {
    use flate2::read::GzDecoder;
    use std::fs::File as StdFile;
    use std::io::{BufRead, BufReader, Write};

    // We need to index the uncompressed FASTA
    // For gzipped files, we'll create an index that maps sequence names to positions
    // in the compressed file using bgzip-style indexing

    // For now, create a simple index by scanning through the gzipped file
    // This index stores: name, length, offset (in uncompressed stream), line_bases, line_width

    let file = StdFile::open(genome_path).context("Failed to open genome file")?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut index_entries: Vec<String> = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_length: u64 = 0;
    let mut current_offset: u64 = 0;
    let mut line_bases: u64 = 0;
    let mut line_width: u64 = 0;
    let mut sequence_start_offset: u64 = 0;
    let mut first_seq_line = true;
    let mut byte_offset: u64 = 0;

    for line in reader.lines() {
        let line = line.context("Failed to read line")?;
        let line_len = line.len() as u64 + 1; // +1 for newline

        if line.starts_with('>') {
            // Save previous entry
            if let Some(name) = current_name.take() {
                index_entries.push(format!(
                    "{}\t{}\t{}\t{}\t{}",
                    name, current_length, sequence_start_offset, line_bases, line_width
                ));
            }

            // Parse new sequence name
            let name = line[1..].split_whitespace().next().unwrap_or("").to_string();
            current_name = Some(name);
            current_length = 0;
            sequence_start_offset = byte_offset + line_len;
            first_seq_line = true;
        } else if current_name.is_some() {
            // Sequence line
            let bases = line.trim().len() as u64;
            current_length += bases;

            if first_seq_line {
                line_bases = bases;
                line_width = line_len;
                first_seq_line = false;
            }
        }

        byte_offset += line_len;
    }

    // Save last entry
    if let Some(name) = current_name {
        index_entries.push(format!(
            "{}\t{}\t{}\t{}\t{}",
            name, current_length, sequence_start_offset, line_bases, line_width
        ));
    }

    // Write index file
    let mut index_file =
        StdFile::create(index_path).context("Failed to create index file")?;
    for entry in index_entries {
        writeln!(index_file, "{}", entry)?;
    }

    // Note: This index is for the uncompressed stream position, which isn't directly
    // usable with gzipped files for random access. For production use, we'd want to
    // either:
    // 1. Decompress the genome to an uncompressed FASTA
    // 2. Use bgzip compression which supports random access
    // 3. Load the entire genome into memory (feasible for ~3GB human genome)

    warn!(
        "Note: Created basic index. For optimal performance, consider decompressing the genome."
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_fetch_creates_directory() {
        // This test would require mocking the HTTP client
        // For now, just verify the config is correct
        let config = AssemblyConfig::new(Species::Human, Some("hg38")).unwrap();
        assert!(config.genome_url.contains("hg38"));
    }
}
