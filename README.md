# Primer Genie

A high-performance Rust CLI tool for designing gene-specific PCR primers with built-in specificity checking.

## Features

- **Automatic Reference Download**: Fetches genome sequences and gene annotations from UCSC
- **Multi-Species Support**: Human (hg38, hg19) and Mouse (mm39, mm10)
- **Inter-Exonic Primer Design**: Primers span intron-exon boundaries to distinguish cDNA from genomic DNA
- **Built-in Specificity Checking**: Filters out primers that match multiple genomic targets
- **Parallel Processing**: Utilizes all CPU cores for fast primer generation
- **Portable Output**: JSONL format for easy downstream processing

## Installation

### From Source

```bash
# Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Clone and build
git clone https://github.com/yourusername/primer-genie.git
cd primer-genie
cargo build --release

# Binary will be at ./target/release/primer-genie
```

## Quick Start

```bash
# 1. Download reference files (one-time setup)
primer-genie fetch human

# 2. Generate primers for all genes
primer-genie generate human -o human_primers.jsonl

# 3. Generate primers for specific genes only
primer-genie generate human --genes my_genes.txt -o subset.jsonl
```

## Usage

### Fetch Reference Files

Downloads genome FASTA and gene annotations from UCSC:

```bash
primer-genie fetch <species> [--assembly <version>]

# Examples:
primer-genie fetch human                    # Downloads hg38 (default)
primer-genie fetch human --assembly hg19    # Downloads hg19
primer-genie fetch mouse                    # Downloads mm39 (default)
primer-genie fetch mouse --assembly mm10    # Downloads mm10
```

Files are cached in `~/.local/share/primer-genie/cache/` (Linux) or equivalent XDG directory.

### Generate Primers

```bash
primer-genie generate <species> [options]

Options:
  -o, --output <file>       Output file (default: stdout)
  -a, --assembly <version>  Genome assembly (default: hg38/mm39)
  --primer-length <n>       Primer length in bp (default: 20)
  --min-tm <n>              Minimum melting temp °C (default: 50.0)
  --max-tm <n>              Maximum melting temp °C (default: 65.0)
  --min-gc <n>              Minimum GC content % (default: 50.0)
  --max-gc <n>              Maximum GC content % (default: 65.0)
  --max-homopolymer <n>     Max consecutive identical bases (default: 3)
  --genes <file>            Only process genes in file (one per line)
  -t, --threads <n>         Worker threads (default: all CPUs)
  --skip-specificity        Skip specificity check (faster)
```

### Examples

```bash
# Full human genome, all genes (~50 min)
primer-genie generate human -o human_primers.jsonl

# Quick test with specificity skipped (~30 min)
primer-genie generate human -o quick.jsonl --skip-specificity

# Custom primer parameters
primer-genie generate human \
  --primer-length 22 \
  --min-tm 55 \
  --max-tm 62 \
  --min-gc 45 \
  --max-gc 60 \
  -o custom_primers.jsonl

# Subset of genes
echo -e "TP53\nBRCA1\nEGFR" > genes.txt
primer-genie generate human --genes genes.txt -o cancer_genes.jsonl
```

## Output Format

Output is JSONL (one JSON object per line):

```json
{"gene":"TP53","transcript":"NM_000546","species":"human","assembly":"hg38","chromosome":"chr17","forward_seq":"GCTATCGGATCCGCAGCCCC","reverse_seq":"TGAAGAACATCGCTGTAGAT","forward_tm":58.42,"reverse_tm":59.18,"forward_gc":65.0,"reverse_gc":55.0,"product_size":85,"introns_spanned":1,"forward_pos":7675000,"reverse_pos":7675085,"targets":1}
```

| Field | Type | Description |
|-------|------|-------------|
| `gene` | string | Gene symbol (e.g., "TP53") |
| `transcript` | string | RefSeq transcript accession |
| `species` | string | Species name |
| `assembly` | string | Genome assembly version |
| `chromosome` | string | Chromosome name |
| `forward_seq` | string | Forward primer sequence (5'→3') |
| `reverse_seq` | string | Reverse primer sequence (5'→3') |
| `forward_tm` | float | Forward primer melting temperature (°C) |
| `reverse_tm` | float | Reverse primer melting temperature (°C) |
| `forward_gc` | float | Forward primer GC content (%) |
| `reverse_gc` | float | Reverse primer GC content (%) |
| `product_size` | int | Expected amplicon size (bp) |
| `introns_spanned` | int | Number of introns between primers |
| `forward_pos` | int | Genomic position of forward primer |
| `reverse_pos` | int | Genomic position of reverse primer |
| `targets` | int | Number of matching genomic locations (1 = specific) |

### Filtering for Specific Primers

For production use, filter to keep only specific primers:

```bash
# Using jq
jq -c 'select(.targets == 1)' primers.jsonl > specific_primers.jsonl

# Using grep (faster for large files)
grep '"targets":1' primers.jsonl > specific_primers.jsonl
```

## Algorithm Details

### Primer Design Strategy

1. **Parse Annotations**: Load RefSeq gene annotations (refGene.txt)
2. **Extract Sequences**: Get spliced transcript sequences (concatenated exons)
3. **Find Candidates**: Scan for primer positions that span exon junctions
4. **Filter Quality**: Apply Tm, GC%, homopolymer, and self-complementarity filters
5. **Check Specificity**: Count how many transcripts each primer pair could amplify

### Thermodynamics

Melting temperature is calculated using the **SantaLucia nearest-neighbor method**:
- Nearest-neighbor enthalpy/entropy parameters (SantaLucia 1998)
- Salt correction for 50mM Na⁺
- Primer concentration of 250nM

### Specificity Checking

Uses **Aho-Corasick** multi-pattern matching to efficiently find all primer matches across the transcriptome. A primer pair is considered to "match" a transcript if:
- Forward primer matches the transcript sequence
- Reverse primer (reverse complement) matches downstream
- Distance between matches is ≤500bp

## Performance

Typical run times on a modern multi-core system:

| Species | Genes | Phase 1 (Design) | Phase 2 (Specificity) | Total |
|---------|-------|------------------|----------------------|-------|
| Human (hg38) | 28,307 | ~30 min | ~23 min | ~53 min |
| Mouse (mm39) | ~25,000 | ~25 min | ~20 min | ~45 min |

Memory usage: ~4-6 GB (genome loaded into memory)

## Reference Data Sources

All reference data is downloaded from UCSC Genome Browser:

| Species | Assembly | Genome | Annotations |
|---------|----------|--------|-------------|
| Human | hg38 | [hg38.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) | [refGene.txt.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz) |
| Human | hg19 | [hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) | [refGene.txt.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz) |
| Mouse | mm39 | [mm39.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz) | [refGene.txt.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refGene.txt.gz) |
| Mouse | mm10 | [mm10.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz) | [refGene.txt.gz](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz) |

## Development

```bash
# Run tests
cargo test

# Run with debug logging
RUST_LOG=debug primer-genie generate human --genes test.txt

# Check code without building
cargo check

# Format code
cargo fmt

# Lint
cargo clippy
```

## License

MIT License - See [LICENSE](LICENSE) for details.

## References

- SantaLucia J Jr. (1998) "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95:1460-1465
- UCSC Genome Browser: https://genome.ucsc.edu/
