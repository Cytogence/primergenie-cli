use primer_genie::primer::thermodynamics::{calculate_gc_content, calculate_tm, longest_homopolymer};

#[test]
fn test_gc_content_balanced() {
    // 50% GC
    let gc = calculate_gc_content(b"ATCGATCGATCGATCGATCG");
    assert!((gc - 50.0).abs() < 0.1, "Expected 50%, got {}", gc);
}

#[test]
fn test_gc_content_gc_rich() {
    // 100% GC
    let gc = calculate_gc_content(b"GCGCGCGCGCGCGCGCGCGC");
    assert!((gc - 100.0).abs() < 0.1, "Expected 100%, got {}", gc);
}

#[test]
fn test_gc_content_at_rich() {
    // 0% GC
    let gc = calculate_gc_content(b"ATATATATATATATATATAT");
    assert!((gc - 0.0).abs() < 0.1, "Expected 0%, got {}", gc);
}

#[test]
fn test_tm_basic() {
    // Standard 20-mer with ~55% GC should have Tm around 55-65°C
    let tm = calculate_tm(b"GCTATCGGATCCGCAGCCCC");
    assert!(tm > 50.0 && tm < 75.0, "Tm {} out of expected range", tm);
}

#[test]
fn test_tm_gc_rich_higher() {
    // GC-rich primers should have higher Tm than AT-rich
    let tm_gc = calculate_tm(b"GCGCGCGCGCGCGCGCGCGC");
    let tm_at = calculate_tm(b"ATATATATATATATATATAT");
    assert!(
        tm_gc > tm_at,
        "GC-rich Tm ({}) should be higher than AT-rich ({})",
        tm_gc,
        tm_at
    );
}

#[test]
fn test_tm_length_effect() {
    // Longer primers should have higher Tm (all else equal)
    let tm_short = calculate_tm(b"ATCGATCGATCG");
    let tm_long = calculate_tm(b"ATCGATCGATCGATCGATCG");
    assert!(
        tm_long > tm_short,
        "Longer primer Tm ({}) should be higher than shorter ({})",
        tm_long,
        tm_short
    );
}

#[test]
fn test_homopolymer_none() {
    let hp = longest_homopolymer(b"ATCGATCGATCG");
    assert_eq!(hp, 1, "No homopolymer, expected 1");
}

#[test]
fn test_homopolymer_run_of_4() {
    let hp = longest_homopolymer(b"ATCGAAAATCG");
    assert_eq!(hp, 4, "Run of 4 A's, expected 4");
}

#[test]
fn test_homopolymer_multiple_runs() {
    let hp = longest_homopolymer(b"AAATTTCCCGGG");
    assert_eq!(hp, 3, "Multiple runs of 3, expected 3");
}

// Test known sequences with expected Tm values
// These can be validated against online Tm calculators or lab data
#[test]
fn test_known_sequences() {
    // These are approximate tests - real validation would compare to
    // experimental or reference Tm values

    // TP53 forward primer (from original database)
    let tp53_fwd = b"GCTATCGGATCCGCAGCCCC";
    let tm = calculate_tm(tp53_fwd);
    // Expected range based on standard calculations
    assert!(
        tm > 55.0 && tm < 70.0,
        "TP53 forward Tm {} outside expected range 55-70",
        tm
    );

    // GC content check
    let gc = calculate_gc_content(tp53_fwd);
    // 13 G/C out of 20 = 65%
    assert!((gc - 65.0).abs() < 0.1, "Expected 65% GC, got {}", gc);
}
