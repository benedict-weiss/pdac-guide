"""Discrimination scoring (primary objective) and an on-target activity proxy."""

from __future__ import annotations

import re
from dataclasses import dataclass

from .cfd import cfd_score
from .enumerate import Guide

# Spacer positions >= SEED_START (1-based, PAM-proximal) count as "seed".
SEED_START = 11


@dataclass(frozen=True)
class DiscriminationResult:
    score_mut: float
    score_wt: float
    margin: float
    mechanism: str   # 'pam_differential' | 'seed_mismatch' | 'distal_mismatch' | 'none'


def discrimination(guide: Guide) -> DiscriminationResult:
    score_mut = 1.0   # gRNA matches the mutant protospacer by construction
    if not guide.wt_has_pam:
        score_wt = 0.0
        mechanism = "pam_differential"
    else:
        score_wt = cfd_score(guide.spacer, guide.wt_protospacer)
        if guide.snv_pos is None:
            mechanism = "none"            # SNV only in PAM but WT also has a PAM
        elif guide.snv_pos >= SEED_START:
            mechanism = "seed_mismatch"
        else:
            mechanism = "distal_mismatch"
    return DiscriminationResult(score_mut, score_wt, score_mut - score_wt, mechanism)


def on_target_proxy(spacer: str) -> float:
    """Simple, documented activity guardrail in [0, 1].

    NOT a validated efficiency predictor. Published predictors (Rule Set 2 etc.)
    are trained on Cas9 *cutting*; this trigger is *binding*. Used only as a
    tie-break among already-safe, well-discriminating guides.
    """
    s = spacer.upper()
    gc = (s.count("G") + s.count("C")) / len(s)
    score = 1.0
    if gc < 0.40 or gc > 0.70:
        score *= 0.7
    if re.search(r"(.)\1{3,}", s):     # any homopolymer run >= 4
        score *= 0.7
    if "TTTT" in s:                    # Pol III terminator
        score *= 0.5
    return score
