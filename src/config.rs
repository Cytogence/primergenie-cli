use anyhow::{anyhow, Result};
use directories::ProjectDirs;
use std::path::PathBuf;

/// Species supported by primer-genie
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Species {
    Human,
    Mouse,
}

impl Species {
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "human" | "homo_sapiens" | "h_sapiens" => Ok(Species::Human),
            "mouse" | "mus_musculus" | "m_musculus" => Ok(Species::Mouse),
            _ => Err(anyhow!("Unknown species: {}. Supported: human, mouse", s)),
        }
    }

    pub fn name(&self) -> &'static str {
        match self {
            Species::Human => "human",
            Species::Mouse => "mouse",
        }
    }

    pub fn default_assembly(&self) -> &'static str {
        match self {
            Species::Human => "hg38",
            Species::Mouse => "mm39",
        }
    }

    pub fn valid_assemblies(&self) -> &'static [&'static str] {
        match self {
            Species::Human => &["hg38", "hg19"],
            Species::Mouse => &["mm39", "mm10"],
        }
    }
}

/// Assembly configuration with download URLs
#[derive(Debug, Clone)]
pub struct AssemblyConfig {
    pub species: Species,
    pub assembly: String,
    pub genome_url: String,
    pub annotations_url: String,
}

impl AssemblyConfig {
    pub fn new(species: Species, assembly: Option<&str>) -> Result<Self> {
        let assembly = assembly.unwrap_or(species.default_assembly());

        // Validate assembly for species
        if !species.valid_assemblies().contains(&assembly) {
            return Err(anyhow!(
                "Invalid assembly '{}' for {}. Valid assemblies: {:?}",
                assembly,
                species.name(),
                species.valid_assemblies()
            ));
        }

        let base_url = format!(
            "https://hgdownload.soe.ucsc.edu/goldenPath/{}/",
            assembly
        );

        let genome_url = format!("{}bigZips/{}.fa.gz", base_url, assembly);
        let annotations_url = format!("{}database/refGene.txt.gz", base_url);

        Ok(Self {
            species,
            assembly: assembly.to_string(),
            genome_url,
            annotations_url,
        })
    }
}

/// Get the data directory for primer-genie
pub fn data_dir() -> Result<PathBuf> {
    if let Some(proj_dirs) = ProjectDirs::from("com", "primer-genie", "primer-genie") {
        Ok(proj_dirs.data_dir().to_path_buf())
    } else {
        // Fallback to ~/.primer-genie
        let home = dirs_home()?;
        Ok(home.join(".primer-genie"))
    }
}

/// Get the cache directory for a specific species/assembly
pub fn cache_dir(species: Species, assembly: &str) -> Result<PathBuf> {
    let base = data_dir()?;
    Ok(base.join("cache").join(species.name()).join(assembly))
}

/// Get the genome file path
pub fn genome_path(species: Species, assembly: &str) -> Result<PathBuf> {
    Ok(cache_dir(species, assembly)?.join("genome.fa.gz"))
}

/// Get the FASTA index file path
pub fn genome_index_path(species: Species, assembly: &str) -> Result<PathBuf> {
    Ok(cache_dir(species, assembly)?.join("genome.fa.gz.fai"))
}

/// Get the gene index file path (bgzip + tabix indexed)
pub fn genome_bgzip_path(species: Species, assembly: &str) -> Result<PathBuf> {
    Ok(cache_dir(species, assembly)?.join("genome.fa.bgz"))
}

/// Get the annotations file path
pub fn annotations_path(species: Species, assembly: &str) -> Result<PathBuf> {
    Ok(cache_dir(species, assembly)?.join("refGene.txt.gz"))
}

/// Get home directory
fn dirs_home() -> Result<PathBuf> {
    dirs::home_dir().ok_or_else(|| anyhow!("Could not determine home directory"))
}

/// Check if reference files exist for a species/assembly
pub fn references_exist(species: Species, assembly: &str) -> Result<bool> {
    let genome = genome_path(species, assembly)?;
    let annotations = annotations_path(species, assembly)?;

    Ok(genome.exists() && annotations.exists())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_species_from_str() {
        assert!(matches!(Species::from_str("human"), Ok(Species::Human)));
        assert!(matches!(Species::from_str("Human"), Ok(Species::Human)));
        assert!(matches!(Species::from_str("mouse"), Ok(Species::Mouse)));
        assert!(Species::from_str("unknown").is_err());
    }

    #[test]
    fn test_assembly_config() {
        let config = AssemblyConfig::new(Species::Human, None).unwrap();
        assert_eq!(config.assembly, "hg38");
        assert!(config.genome_url.contains("hg38"));

        let config = AssemblyConfig::new(Species::Human, Some("hg19")).unwrap();
        assert_eq!(config.assembly, "hg19");

        assert!(AssemblyConfig::new(Species::Human, Some("mm39")).is_err());
    }
}
