from guide_design.targets import Target
from guide_design.enumerate import enumerate_guides
from guide_design.cfd import cfd_score, POSITION_WEIGHTS
from guide_design.score import discrimination, on_target_proxy


def test_cfd_perfect_match_is_one():
    assert cfd_score("A" * 20, "A" * 20) == 1.0


def test_cfd_seed_mismatch_lower_than_distal():
    spacer = "A" * 20
    distal = "C" + "A" * 19    # mismatch at index 0 (PAM-distal)
    seed = "A" * 19 + "C"      # mismatch at index 19 (PAM-proximal)
    assert cfd_score(spacer, distal) > cfd_score(spacer, seed)


def test_position_weights_monotone_nonincreasing():
    assert len(POSITION_WEIGHTS) == 20
    assert all(POSITION_WEIGHTS[i] >= POSITION_WEIGHTS[i + 1] for i in range(19))


def _target(ref, idx, rb, ab, allele="X"):
    return Target("T", "G", allele, ref, idx, rb, ab, 1.0)


def test_discrimination_seed_mechanism_high_margin():
    # SNV at PAM-proximal end -> seed mismatch -> high margin.
    t = _target("A" * 20 + "TGGAAA", 19, "A", "G")
    g = next(x for x in enumerate_guides(t) if x.strand == "+" and x.snv_pos == 20)
    d = discrimination(g)
    assert d.score_mut == 1.0
    assert d.margin > 0.9
    assert d.mechanism == "seed_mismatch"


def test_discrimination_mutant_specific_pam_zero_wt():
    t = _target("A" * 20 + "TAGAAA", 21, "A", "G")
    g = next(x for x in enumerate_guides(t) if x.mutant_specific_pam)
    d = discrimination(g)
    assert d.score_wt == 0.0
    assert d.margin == 1.0
    assert d.mechanism == "pam_differential"


def test_on_target_proxy_penalizes_polyt():
    assert on_target_proxy("AAGCTAGCTAGCTAGCTAGC") > on_target_proxy("TTTTAGCTAGCTAGCTAGCT")
