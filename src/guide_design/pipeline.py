"""End-to-end orchestration: targets -> guides -> scores -> recommendations."""

from __future__ import annotations

import pandas as pd

from .coverage import coverage_curve
from .enumerate import enumerate_guides
from .rank import RankedGuide, evaluate_guide, recommend_per_allele
from .targets import Target


def run(
    targets: list[Target],
    reference: str | None,
    margin_threshold: float = 0.5,
) -> tuple[dict[str, RankedGuide | None], pd.DataFrame, list[dict]]:
    rgs: list[RankedGuide] = []
    for t in targets:
        for g in enumerate_guides(t):
            rgs.append(evaluate_guide(g, reference=reference))

    recommendations = recommend_per_allele(rgs, margin_threshold)

    rows = []
    for t in targets:
        rec = recommendations.get(t.allele)
        if rec is None:
            rows.append({
                "allele": t.allele, "undesignable": True,
                "recommended_spacer": None, "strand": None, "pam": None,
                "mechanism": None, "margin": None, "score_wt": None,
                "specificity": None, "on_target": None, "frequency": t.frequency,
            })
        else:
            rows.append({
                "allele": t.allele, "undesignable": False,
                "recommended_spacer": rec.guide.spacer, "strand": rec.guide.strand,
                "pam": rec.guide.pam, "mechanism": rec.mechanism,
                "margin": round(rec.margin, 4), "score_wt": round(rec.score_wt, 4),
                "specificity": None if rec.specificity is None else round(rec.specificity, 2),
                "on_target": round(rec.on_target, 4), "frequency": t.frequency,
            })
    table = pd.DataFrame(rows).sort_values("frequency", ascending=False).reset_index(drop=True)
    cov = coverage_curve(recommendations, targets)
    return recommendations, table, cov
