from guide_design.offtarget import (
    OffTarget, find_offtargets, specificity_score, passes_safety_veto,
)


def test_finds_exact_on_target_site():
    spacer = "ACGTACGTACGTACGTACGT"
    reference = "TTT" + spacer + "TGG" + "TTT"   # spacer + NGG PAM
    hits = find_offtargets(spacer, reference, max_mismatch=4)
    exact = [h for h in hits if h.mismatches == 0]
    assert len(exact) == 1
    assert exact[0].strand == "+"


def test_finds_near_match_within_threshold_only():
    spacer = "ACGTACGTACGTACGTACGT"
    near = "ACGTACGTACGTACGTACGA"          # 1 mismatch at PAM-proximal end
    far = "TTTTTTTTTTTTTTTTTTTT"           # >4 mismatches
    reference = "AA" + near + "AGG" + "AA" + far + "CGG"
    hits = find_offtargets(spacer, reference, max_mismatch=4)
    seqs = {h.site_seq for h in hits}
    assert near in seqs
    assert far not in seqs


def test_specificity_score_decreases_with_offtargets():
    assert specificity_score([]) == 100.0
    assert specificity_score([0.1]) < 100.0
    assert specificity_score([0.8]) < specificity_score([0.1])


def test_safety_veto_fails_on_high_cfd_offtarget():
    good = [OffTarget("+", 10, "A" * 20, 3, 0.2)]
    bad = [OffTarget("+", 10, "A" * 20, 1, 0.8)]
    assert passes_safety_veto(good, cfd_threshold=0.5) is True
    assert passes_safety_veto(bad, cfd_threshold=0.5) is False
