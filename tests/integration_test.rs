use primer_genie::config::{AssemblyConfig, Species};
use primer_genie::genome::transcript::{Exon, Strand, Transcript};
use primer_genie::primer::{PrimerParams, PrimerFilter};
use primer_genie::output::PrimerOutput;

#[test]
fn test_species_parsing() {
    assert!(matches!(Species::from_str("human"), Ok(Species::Human)));
    assert!(matches!(Species::from_str("Human"), Ok(Species::Human)));
    assert!(matches!(Species::from_str("HUMAN"), Ok(Species::Human)));
    assert!(matches!(Species::from_str("mouse"), Ok(Species::Mouse)));
    assert!(matches!(Species::from_str("Mouse"), Ok(Species::Mouse)));
    assert!(Species::from_str("unknown").is_err());
}

#[test]
fn test_assembly_config_human() {
    let config = AssemblyConfig::new(Species::Human, None).unwrap();
    assert_eq!(config.assembly, "hg38");
    assert!(config.genome_url.contains("hg38"));
    assert!(config.annotations_url.contains("refGene.txt.gz"));

    let config = AssemblyConfig::new(Species::Human, Some("hg19")).unwrap();
    assert_eq!(config.assembly, "hg19");
    assert!(config.genome_url.contains("hg19"));

    // Invalid assembly for human
    assert!(AssemblyConfig::new(Species::Human, Some("mm39")).is_err());
}

#[test]
fn test_assembly_config_mouse() {
    let config = AssemblyConfig::new(Species::Mouse, None).unwrap();
    assert_eq!(config.assembly, "mm39");
    assert!(config.genome_url.contains("mm39"));

    let config = AssemblyConfig::new(Species::Mouse, Some("mm10")).unwrap();
    assert_eq!(config.assembly, "mm10");

    // Invalid assembly for mouse
    assert!(AssemblyConfig::new(Species::Mouse, Some("hg38")).is_err());
}

#[test]
fn test_transcript_creation() {
    let transcript = Transcript {
        accession: "NM_000546".to_string(),
        gene_symbol: "TP53".to_string(),
        chromosome: "chr17".to_string(),
        strand: Strand::Minus,
        tx_start: 7668402,
        tx_end: 7687550,
        cds_start: 7669609,
        cds_end: 7676594,
        exons: vec![
            Exon::new(7668402, 7669690, 11),
            Exon::new(7670609, 7670715, 10),
            Exon::new(7673534, 7673608, 9),
        ],
    };

    assert_eq!(transcript.accession, "NM_000546");
    assert_eq!(transcript.gene_symbol, "TP53");
    assert!(transcript.is_coding());
    assert_eq!(transcript.intron_count(), 2);
}

#[test]
fn test_primer_filter() {
    let params = PrimerParams {
        primer_length: 20,
        min_tm: 50.0,
        max_tm: 70.0,
        min_gc: 40.0,
        max_gc: 70.0,
        max_homopolymer: 4,
        ..PrimerParams::default()
    };
    let filter = PrimerFilter::new(params);

    // Good primer
    assert!(filter.passes(b"GCTATCGGATCCGCAGCCCC"));

    // Contains N - should fail
    assert!(!filter.passes(b"GCTATCGNNTCCGCAGCCCC"));

    // Wrong length - should fail
    assert!(!filter.passes(b"ATCG"));

    // Homopolymer run of 5 - should fail
    assert!(!filter.passes(b"ATCGAAAAATCGCGCGCGCG"));
}

#[test]
fn test_primer_output_serialization() {
    let output = PrimerOutput {
        gene: "TP53".to_string(),
        transcript: "NM_000546".to_string(),
        species: "human".to_string(),
        assembly: "hg38".to_string(),
        chromosome: "chr17".to_string(),
        forward_seq: "GCTATCGGATCCGCAGCCCC".to_string(),
        reverse_seq: "TGAAGAACATCGCTGTAGAT".to_string(),
        forward_tm: 58.42,
        reverse_tm: 59.18,
        forward_gc: 65.0,
        reverse_gc: 55.0,
        product_size: 85,
        introns_spanned: 1,
        forward_pos: 7675000,
        reverse_pos: 7675085,
        targets: 1,
    };

    let json = serde_json::to_string(&output).unwrap();

    // Verify JSON contains expected fields
    assert!(json.contains("\"gene\":\"TP53\""));
    assert!(json.contains("\"species\":\"human\""));
    assert!(json.contains("\"assembly\":\"hg38\""));
    assert!(json.contains("\"targets\":1"));
    assert!(json.contains("\"product_size\":85"));
    assert!(json.contains("\"introns_spanned\":1"));

    // Verify deserialization
    let parsed: PrimerOutput = serde_json::from_str(&json).unwrap();
    assert_eq!(parsed.gene, "TP53");
    assert_eq!(parsed.targets, 1);
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
fn test_exon_length() {
    let exon = Exon::new(1000, 1500, 1);
    assert_eq!(exon.length(), 500);
}
