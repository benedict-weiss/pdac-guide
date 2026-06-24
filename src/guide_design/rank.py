"""Combine discrimination, off-target safety, and on-target proxy into a ranking."""

from __future__ import annotations

from dataclasses import dataclass

from .enumerate import Guide
from .offtarget import find_offtargets, passes_safety_veto, specificity_score
from .score import discrimination, on_target_proxy


@dataclass(frozen=True)
class RankedGuide:
    guide: Guide
    margin: float
    score_wt: float
    mechanism: str
    specificity: float | None
    on_target: float
    veto_pass: bool


def evaluate_guide(
    guide: Guide,
    reference: str | None = None,
    max_mismatch: int = 4,
    cfd_threshold: float = 0.5,
) -> RankedGuide:
    d = discrimination(guide)
    if reference is None:
        specificity: float | None = None
        veto = True
    else:
        offs = [o for o in find_offtargets(guide.spacer, reference, max_mismatch=max_mismatch)
                if o.mismatches > 0]
        specificity = specificity_score([o.cfd for o in offs])
        veto = passes_safety_veto(offs, cfd_threshold)
    return RankedGuide(
        guide=guide,
        margin=d.margin,
        score_wt=d.score_wt,
        mechanism=d.mechanism,
        specificity=specificity,
        on_target=on_target_proxy(guide.spacer),
        veto_pass=veto,
    )


def _qualifies(rg: RankedGuide, margin_threshold: float) -> bool:
    return rg.veto_pass and rg.margin >= margin_threshold


def rank_guides(rgs: list[RankedGuide], margin_threshold: float = 0.5) -> list[RankedGuide]:
    return sorted(
        rgs,
        key=lambda r: (_qualifies(r, margin_threshold), r.margin, r.on_target),
        reverse=True,
    )


def recommend_per_allele(
    rgs: list[RankedGuide], margin_threshold: float = 0.5
) -> dict[str, RankedGuide | None]:
    out: dict[str, RankedGuide | None] = {}
    for allele in sorted({r.guide.allele for r in rgs}):
        ranked = rank_guides([r for r in rgs if r.guide.allele == allele], margin_threshold)
        out[allele] = ranked[0] if ranked and _qualifies(ranked[0], margin_threshold) else None
    return out
