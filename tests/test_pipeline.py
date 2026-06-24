import pandas as pd

from guide_design.targets import Target
from guide_design.pipeline import run


def _designable(allele, freq):
    return Target("KRAS_" + allele, "KRAS", allele, "A" * 20 + "TGGAAA", 19, "A", "G", freq)


def test_pipeline_returns_row_per_allele_and_recommends_designable():
    targets = [_designable("G12D", 0.40), _designable("G12V", 0.30)]
    # Small synthetic reference so the off-target path runs deterministically.
    reference = "GGGTTT" + "A" * 20 + "TGG" + "CCCCCC"
    rec, table, cov = run(targets, reference=reference, margin_threshold=0.5)
    assert isinstance(table, pd.DataFrame)
    assert set(table["allele"]) == {"G12D", "G12V"}
    assert (~table["undesignable"]).all()
    assert table["recommended_spacer"].notna().all()
    assert len(cov) == 2


def test_pipeline_flags_undesignable_allele():
    # Distal-only mismatch + high threshold -> undesignable.
    t = Target("KRAS_G12X", "KRAS", "G12X", "C" + "A" * 19 + "TGGAAA", 0, "C", "G", 0.1)
    rec, table, cov = run([t], reference=None, margin_threshold=0.9)
    assert bool(table.loc[table["allele"] == "G12X", "undesignable"].iloc[0]) is True
    assert cov == []
