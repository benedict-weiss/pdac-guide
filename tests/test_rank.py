from guide_design.targets import Target
from guide_design.enumerate import enumerate_guides
from guide_design.rank import evaluate_guide, rank_guides, recommend_per_allele


def _target(ref, idx, rb, ab, allele="X", freq=1.0):
    return Target("KRAS_" + allele, "KRAS", allele, ref, idx, rb, ab, freq)


def test_evaluate_guide_without_reference_skips_offtarget():
    t = _target("A" * 20 + "TGGAAA", 19, "A", "G", allele="G12D")
    g = enumerate_guides(t)[0]
    rg = evaluate_guide(g, reference=None)
    assert rg.specificity is None
    assert rg.veto_pass is True
    assert rg.margin > 0.9


def test_rank_orders_qualified_first_then_by_margin():
    t_seed = _target("A" * 20 + "TGGAAA", 19, "A", "G", allele="A")   # seed -> high margin
    t_distal = _target("C" + "A" * 19 + "TGGAAA", 0, "C", "G", allele="B")  # distal -> low margin
    rgs = [evaluate_guide(enumerate_guides(t_distal)[0]),
           evaluate_guide(enumerate_guides(t_seed)[0])]
    ranked = rank_guides(rgs, margin_threshold=0.5)
    assert ranked[0].margin > ranked[1].margin    # high-margin seed guide first


def test_recommend_marks_undesignable_allele_none():
    # Only a distal mismatch available -> margin below threshold -> undesignable.
    t = _target("C" + "A" * 19 + "TGGAAA", 0, "C", "G", allele="G12X")
    rgs = [evaluate_guide(g) for g in enumerate_guides(t)]
    rec = recommend_per_allele(rgs, margin_threshold=0.9)
    assert rec["G12X"] is None
