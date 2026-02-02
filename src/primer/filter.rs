use super::thermodynamics::{
    calculate_gc_content, calculate_tm, longest_homopolymer, self_complementarity,
};

/// Parameters for primer filtering
#[derive(Debug, Clone)]
pub struct PrimerParams {
    /// Primer length in base pairs
    pub primer_length: usize,
    /// Minimum melting temperature (°C)
    pub min_tm: f64,
    /// Maximum melting temperature (°C)
    pub max_tm: f64,
    /// Minimum GC content (%)
    pub min_gc: f64,
    /// Maximum GC content (%)
    pub max_gc: f64,
    /// Maximum consecutive identical bases
    pub max_homopolymer: usize,
    /// Maximum self-complementarity length
    pub max_self_complementarity: usize,
    /// Maximum 3' complementarity (for primer dimers)
    pub max_3prime_complementarity: usize,
}

impl Default for PrimerParams {
    fn default() -> Self {
        Self {
            primer_length: 20,
            min_tm: 50.0,
            max_tm: 65.0,
            min_gc: 50.0,
            max_gc: 65.0,
            max_homopolymer: 3,
            max_self_complementarity: 4,
            max_3prime_complementarity: 3,
        }
    }
}

/// Result of primer filtering with reason for rejection
#[derive(Debug, Clone)]
pub enum FilterResult {
    Pass,
    FailTmLow(f64),
    FailTmHigh(f64),
    FailGcLow(f64),
    FailGcHigh(f64),
    FailHomopolymer(usize),
    FailSelfComplementarity(usize),
    FailContainsN,
    FailLength(usize),
}

impl FilterResult {
    pub fn is_pass(&self) -> bool {
        matches!(self, FilterResult::Pass)
    }
}

/// Primer filter that checks sequences against quality criteria
pub struct PrimerFilter {
    params: PrimerParams,
}

impl PrimerFilter {
    pub fn new(params: PrimerParams) -> Self {
        Self { params }
    }

    /// Check if a primer sequence passes all filters
    pub fn check(&self, seq: &[u8]) -> FilterResult {
        // Check for N bases
        if seq.iter().any(|&b| b == b'N' || b == b'n') {
            return FilterResult::FailContainsN;
        }

        // Check length
        if seq.len() != self.params.primer_length {
            return FilterResult::FailLength(seq.len());
        }

        // Check Tm
        let tm = calculate_tm(seq);
        if tm < self.params.min_tm {
            return FilterResult::FailTmLow(tm);
        }
        if tm > self.params.max_tm {
            return FilterResult::FailTmHigh(tm);
        }

        // Check GC content
        let gc = calculate_gc_content(seq);
        if gc < self.params.min_gc {
            return FilterResult::FailGcLow(gc);
        }
        if gc > self.params.max_gc {
            return FilterResult::FailGcHigh(gc);
        }

        // Check homopolymer
        let homopolymer = longest_homopolymer(seq);
        if homopolymer > self.params.max_homopolymer {
            return FilterResult::FailHomopolymer(homopolymer);
        }

        // Check self-complementarity
        let self_comp = self_complementarity(seq);
        if self_comp > self.params.max_self_complementarity {
            return FilterResult::FailSelfComplementarity(self_comp);
        }

        FilterResult::Pass
    }

    /// Check if a primer sequence passes all filters (boolean)
    pub fn passes(&self, seq: &[u8]) -> bool {
        self.check(seq).is_pass()
    }

    /// Check 3' complementarity between forward and reverse primers
    pub fn check_3prime_complementarity(&self, forward: &[u8], reverse: &[u8]) -> bool {
        // Check last N bases of each primer for complementarity
        let check_len = self.params.max_3prime_complementarity + 1;
        let check_len = check_len.min(forward.len()).min(reverse.len());

        let fwd_3prime = &forward[forward.len() - check_len..];
        let rev_3prime = &reverse[reverse.len() - check_len..];

        // Reverse complement of reverse primer 3' end
        let rev_3prime_rc: Vec<u8> = rev_3prime
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'G' | b'g' => b'C',
                b'C' | b'c' => b'G',
                _ => b'N',
            })
            .collect();

        // Count complementary bases
        let mut comp_count = 0;
        for (b1, b2) in fwd_3prime.iter().zip(rev_3prime_rc.iter()) {
            if b1.to_ascii_uppercase() == b2.to_ascii_uppercase() {
                comp_count += 1;
            }
        }

        comp_count <= self.params.max_3prime_complementarity
    }

    /// Get the parameters
    pub fn params(&self) -> &PrimerParams {
        &self.params
    }
}

/// Check if a sequence position is valid for primer placement
///
/// Primers should avoid:
/// - Repetitive regions (low complexity)
/// - Positions with many N bases nearby
pub fn is_valid_primer_position(seq: &[u8], pos: usize, primer_len: usize) -> bool {
    if pos + primer_len > seq.len() {
        return false;
    }

    let primer_seq = &seq[pos..pos + primer_len];

    // Check for N bases
    if primer_seq.iter().any(|&b| b == b'N' || b == b'n') {
        return false;
    }

    true
}

/// Calculate the GC clamp score (G/C at 3' end)
///
/// A good primer should have 1-2 G/C bases in the last 5 bases
/// Returns the count of G/C in the last 5 bases
pub fn gc_clamp(seq: &[u8]) -> usize {
    let last_5 = if seq.len() >= 5 {
        &seq[seq.len() - 5..]
    } else {
        seq
    };

    last_5
        .iter()
        .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_pass() {
        let params = PrimerParams {
            primer_length: 20,
            min_tm: 50.0,
            max_tm: 70.0,
            min_gc: 40.0,
            max_gc: 70.0,
            max_homopolymer: 4,
            max_self_complementarity: 5,
            max_3prime_complementarity: 3,
        };
        let filter = PrimerFilter::new(params);

        // Good primer
        let result = filter.check(b"GCTATCGGATCCGCAGCCCC");
        assert!(result.is_pass(), "Should pass: {:?}", result);
    }

    #[test]
    fn test_filter_fail_n() {
        let params = PrimerParams::default();
        let filter = PrimerFilter::new(params);

        let result = filter.check(b"GCTATCGNNTCCGCAGCCCC");
        assert!(matches!(result, FilterResult::FailContainsN));
    }

    #[test]
    fn test_filter_fail_homopolymer() {
        let params = PrimerParams {
            max_homopolymer: 3,
            ..PrimerParams::default()
        };
        let filter = PrimerFilter::new(params);

        let result = filter.check(b"GCTAAAATCCGCAGCCCCGC");
        assert!(matches!(result, FilterResult::FailHomopolymer(_)));
    }

    #[test]
    fn test_gc_clamp() {
        assert_eq!(gc_clamp(b"ATCGATCGATCGATCGATCG"), 2); // ends in TCG
        assert_eq!(gc_clamp(b"ATCGATCGATCGATCGATAT"), 0); // ends in TATAT
        assert_eq!(gc_clamp(b"ATCGATCGATCGATCGCGCG"), 5); // ends in GCGCG
    }

    #[test]
    fn test_3prime_complementarity() {
        let params = PrimerParams {
            max_3prime_complementarity: 3,
            ..PrimerParams::default()
        };
        let filter = PrimerFilter::new(params);

        // Non-complementary 3' ends
        assert!(filter.check_3prime_complementarity(b"ATCGATCGATCGATCGATCG", b"ATCGATCGATCGATCGATCG"));

        // Highly complementary 3' ends (ATCG and CGAT are reverse complements)
        // This should fail if we check strictly
    }
}
