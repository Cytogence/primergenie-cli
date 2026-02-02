/// Thermodynamic calculations for primer design
///
/// Uses the SantaLucia nearest-neighbor method for Tm calculation
/// Reference: SantaLucia J Jr. (1998) PNAS 95:1460-1465

/// Nearest-neighbor enthalpy values (kcal/mol)
/// Format: 5'-XY-3' / 3'-X'Y'-5'
const NN_ENTHALPY: [[f64; 4]; 4] = [
    // A      C      G      T     (second base)
    [-7.9, -8.4, -7.8, -7.2], // A (first base) - AA/TT, AC/TG, AG/TC, AT/TA
    [-8.5, -8.0, -10.6, -7.8], // C - CA/GT, CC/GG, CG/GC, CT/GA
    [-8.2, -9.8, -8.0, -8.4], // G - GA/CT, GC/CG, GG/CC, GT/CA
    [-7.2, -8.2, -8.5, -7.9], // T - TA/AT, TC/AG, TG/AC, TT/AA
];

/// Nearest-neighbor entropy values (cal/mol·K)
const NN_ENTROPY: [[f64; 4]; 4] = [
    // A       C       G       T
    [-22.2, -22.4, -21.0, -20.4], // A
    [-22.7, -19.9, -27.2, -21.0], // C
    [-22.2, -24.4, -19.9, -22.4], // G
    [-21.3, -22.2, -22.7, -22.2], // T
];

/// Initiation parameters
const INIT_ENTHALPY_GC: f64 = 0.1;  // G/C terminal
const INIT_ENTROPY_GC: f64 = -2.8;
const INIT_ENTHALPY_AT: f64 = 2.3;  // A/T terminal
const INIT_ENTROPY_AT: f64 = 4.1;

/// Convert base to index (A=0, C=1, G=2, T=3)
fn base_to_index(base: u8) -> Option<usize> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Calculate melting temperature using the SantaLucia nearest-neighbor method
///
/// Parameters:
/// - `seq`: Primer sequence (5' to 3')
/// - `na_conc`: Sodium concentration in mM (default: 50)
/// - `primer_conc`: Primer concentration in nM (default: 250)
///
/// Returns: Tm in degrees Celsius
pub fn calculate_tm(seq: &[u8]) -> f64 {
    calculate_tm_with_params(seq, 50.0, 250.0)
}

/// Calculate Tm with custom salt and primer concentrations
pub fn calculate_tm_with_params(seq: &[u8], na_conc: f64, primer_conc: f64) -> f64 {
    if seq.len() < 2 {
        return 0.0;
    }

    let mut enthalpy = 0.0;
    let mut entropy = 0.0;

    // Sum nearest-neighbor contributions
    for i in 0..seq.len() - 1 {
        let (b1, b2) = (seq[i], seq[i + 1]);

        if let (Some(i1), Some(i2)) = (base_to_index(b1), base_to_index(b2)) {
            enthalpy += NN_ENTHALPY[i1][i2];
            entropy += NN_ENTROPY[i1][i2];
        }
    }

    // Add initiation parameters based on terminal bases
    // 5' terminal
    if let Some(i) = base_to_index(seq[0]) {
        if i == 0 || i == 3 {
            // A or T
            enthalpy += INIT_ENTHALPY_AT;
            entropy += INIT_ENTROPY_AT;
        } else {
            enthalpy += INIT_ENTHALPY_GC;
            entropy += INIT_ENTROPY_GC;
        }
    }

    // 3' terminal
    if let Some(i) = base_to_index(seq[seq.len() - 1]) {
        if i == 0 || i == 3 {
            // A or T
            enthalpy += INIT_ENTHALPY_AT;
            entropy += INIT_ENTROPY_AT;
        } else {
            enthalpy += INIT_ENTHALPY_GC;
            entropy += INIT_ENTROPY_GC;
        }
    }

    // Convert entropy from cal to kcal
    entropy /= 1000.0;

    // Salt correction (Owczarzy et al., 2004)
    let na_m = na_conc / 1000.0; // Convert mM to M
    let salt_correction = 0.368 * (seq.len() - 1) as f64 * na_m.ln();
    entropy += salt_correction / 1000.0;

    // Calculate Tm
    // Tm = ΔH / (ΔS + R·ln(Ct/4))
    // R = 1.987 cal/(mol·K)
    let r = 1.987 / 1000.0; // kcal/(mol·K)
    let ct = primer_conc / 1e9; // Convert nM to M

    let tm_kelvin = (enthalpy * 1000.0) / (entropy * 1000.0 + r * 1000.0 * (ct / 4.0).ln());
    let tm_celsius = tm_kelvin - 273.15;

    tm_celsius
}

/// Calculate GC content as a percentage
pub fn calculate_gc_content(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }

    let gc_count = seq
        .iter()
        .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
        .count();

    (gc_count as f64 / seq.len() as f64) * 100.0
}

/// Check for homopolymer runs (consecutive identical bases)
///
/// Returns the length of the longest homopolymer run
pub fn longest_homopolymer(seq: &[u8]) -> usize {
    if seq.is_empty() {
        return 0;
    }

    let mut max_run = 1;
    let mut current_run = 1;
    let mut prev_base = seq[0].to_ascii_uppercase();

    for &base in seq.iter().skip(1) {
        let base = base.to_ascii_uppercase();
        if base == prev_base {
            current_run += 1;
            max_run = max_run.max(current_run);
        } else {
            current_run = 1;
            prev_base = base;
        }
    }

    max_run
}

/// Check for dinucleotide repeats
///
/// Returns the number of consecutive dinucleotide repeats
pub fn longest_dinucleotide_repeat(seq: &[u8]) -> usize {
    if seq.len() < 4 {
        return 0;
    }

    let mut max_repeat = 0;

    for i in 0..seq.len() - 3 {
        let di1 = &seq[i..i + 2];
        let mut repeat_count = 1;

        let mut j = i + 2;
        while j + 1 < seq.len() {
            if seq[j].to_ascii_uppercase() == di1[0].to_ascii_uppercase()
                && seq[j + 1].to_ascii_uppercase() == di1[1].to_ascii_uppercase()
            {
                repeat_count += 1;
                j += 2;
            } else {
                break;
            }
        }

        max_repeat = max_repeat.max(repeat_count);
    }

    max_repeat
}

/// Check for self-complementarity (hairpin potential)
///
/// Returns the length of the longest self-complementary region
pub fn self_complementarity(seq: &[u8]) -> usize {
    let seq_upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    let mut max_comp = 0;

    // Check for hairpin-forming sequences
    for gap in 3..=6 {
        // Loop size
        for i in 0..seq.len().saturating_sub(gap + 2) {
            let mut comp_len = 0;
            let mut j = i;
            let mut k = i + gap + 1;

            while k < seq.len() {
                let b1 = seq_upper[j];
                let b2 = seq_upper[k];

                let complement = match b1 {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'G' => b'C',
                    b'C' => b'G',
                    _ => b'N',
                };

                if b2 == complement {
                    comp_len += 1;
                    if j == 0 {
                        break;
                    }
                    j -= 1;
                    k += 1;
                } else {
                    break;
                }
            }

            max_comp = max_comp.max(comp_len);
        }
    }

    max_comp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_content() {
        assert!((calculate_gc_content(b"ATCG") - 50.0).abs() < 0.01);
        assert!((calculate_gc_content(b"GGCC") - 100.0).abs() < 0.01);
        assert!((calculate_gc_content(b"AATT") - 0.0).abs() < 0.01);
        assert!((calculate_gc_content(b"GCGCGCGCGC") - 100.0).abs() < 0.01);
    }

    #[test]
    fn test_homopolymer() {
        assert_eq!(longest_homopolymer(b"ATCG"), 1);
        assert_eq!(longest_homopolymer(b"AAAA"), 4);
        assert_eq!(longest_homopolymer(b"AATTCC"), 2);
        assert_eq!(longest_homopolymer(b"AAATTTCCC"), 3);
    }

    #[test]
    fn test_tm_calculation() {
        // Test with a known sequence
        // These values should be approximately correct for standard conditions
        let tm = calculate_tm(b"GCTATCGGATCCGCAGCCCC");
        assert!(tm > 50.0 && tm < 75.0, "Tm should be reasonable: {}", tm);

        // Longer primer should have higher Tm
        let tm_short = calculate_tm(b"ATCGATCG");
        let tm_long = calculate_tm(b"ATCGATCGATCGATCG");
        assert!(tm_long > tm_short);

        // GC-rich primer should have higher Tm than AT-rich
        let tm_gc = calculate_tm(b"GCGCGCGCGCGCGCGCGCGC");
        let tm_at = calculate_tm(b"ATATATATATATATATATAT");
        assert!(tm_gc > tm_at);
    }

    #[test]
    fn test_dinucleotide_repeat() {
        assert_eq!(longest_dinucleotide_repeat(b"ATCG"), 1);
        assert_eq!(longest_dinucleotide_repeat(b"ATATATATAT"), 5);
        assert_eq!(longest_dinucleotide_repeat(b"GCGCGC"), 3);
    }
}
