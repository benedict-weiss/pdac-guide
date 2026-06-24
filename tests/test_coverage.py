from guide_design.targets import Target
from guide_design.enumerate import enumerate_guides
from guide_design.rank import evaluate_guide, recommend_per_allele
from guide_design.coverage import coverage_curve


def _designable_target(allele, freq):
    # Seed-position SNV -> a qualifying guide always exists.
    return Target("KRAS_" + allele, "KRAS", allele, "A" * 20 + "TGGAAA", 19, "A", "G", freq)


def test_coverage_accumulates_in_frequency_order():
    targets = [_designable_target("G12D", 0.40),
               _designable_target("G12V", 0.30),
               _designable_target("G12R", 0.17)]
    rgs = [evaluate_guide(g) for t in targets for g in enumerate_guides(t)]
    rec = recommend_per_allele(rgs, margin_threshold=0.5)
    rows = coverage_curve(rec, targets, kras_fraction=0.90)
    assert [r["added_allele"] for r in rows] == ["G12D", "G12V", "G12R"]
    assert abs(rows[-1]["cumulative_kras_frac"] - 0.87) < 1e-9
    assert abs(rows[-1]["cumulative_pdac_frac"] - 0.87 * 0.90) < 1e-9
    assert abs(rows[0]["marginal_pdac"] - 0.40 * 0.90) < 1e-9


def test_coverage_skips_undesignable_alleles():
    targets = [_designable_target("G12D", 0.40)]
    rec = {"G12D": None}
    assert coverage_curve(rec, targets) == []
