"""Translate per-allele recommendations into patient-coverage figures."""

from __future__ import annotations

from .rank import RankedGuide
from .targets import Target


def coverage_curve(
    recommendations: dict[str, RankedGuide | None],
    targets: list[Target],
    kras_fraction: float = 0.90,
) -> list[dict]:
    """Cumulative coverage as guides are added in descending allele frequency.

    kras_fraction = fraction of PDAC tumors that are KRAS-mutant (~0.90).
    """
    freq = {t.allele: t.frequency for t in targets}
    designable = [a for a, r in recommendations.items() if r is not None]
    ordered = sorted(designable, key=lambda a: freq[a], reverse=True)
    rows: list[dict] = []
    cum = 0.0
    for a in ordered:
        cum += freq[a]
        rows.append({
            "added_allele": a,
            "cumulative_kras_frac": cum,
            "cumulative_pdac_frac": cum * kras_fraction,
            "marginal_pdac": freq[a] * kras_fraction,
        })
    return rows
