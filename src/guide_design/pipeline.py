"""End-to-end orchestration: targets -> guides -> scores -> recommendations."""

from __future__ import annotations

from typing import Callable

import pandas as pd

from .coverage import coverage_curve
from .enumerate import enumerate_guides
from .rank import RankedGuide, evaluate_guide, recommend_per_allele
from .targets import Target


def run(
    targets: list[Target],
    reference: str | None,
    margin_threshold: float = 0.5,
    on_progress: Callable[[int, int], None] | None = None,
) -> tuple[dict[str, RankedGuide | None], pd.DataFrame, list[dict]]:
    # Materialize candidates up front so progress has a known denominator. The
    # per-guide off-target scan is the expensive step, so one tick per guide is
    # the natural unit of work.
    candidates = [g for t in targets for g in enumerate_guides(t)]
    total = len(candidates)
    if on_progress:
        on_progress(0, total)
    rgs: list[RankedGuide] = []
    for done, g in enumerate(candidates, start=1):
        rgs.append(evaluate_guide(g, reference=reference))
        if on_progress:
            on_progress(done, total)

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
