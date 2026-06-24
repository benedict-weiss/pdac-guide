"""CFD-style mismatch scoring.

POSITION_WEIGHTS is a documented monotonic approximation of the Doench et al.
2016 Cutting Frequency Determination averages: a mismatch at a PAM-distal
position (index 0) retains most activity; a mismatch in the PAM-proximal seed
(index 19) abolishes it. This is an approximation, not the full per-mismatch
empirical table; see the writeup's limitations.
"""

from __future__ import annotations

# index 0 = PAM-distal (5'), index 19 = PAM-proximal (adjacent to PAM / seed).
POSITION_WEIGHTS: list[float] = [
    0.85, 0.80, 0.78, 0.75, 0.72, 0.70, 0.65, 0.60, 0.55, 0.50,
    0.45, 0.40, 0.32, 0.25, 0.18, 0.12, 0.08, 0.05, 0.03, 0.02,
]


def cfd_score(spacer: str, target20: str) -> float:
    """Product of position weights over mismatched positions; 1.0 if identical."""
    if len(spacer) != len(target20):
        raise ValueError("spacer and target must be the same length")
    score = 1.0
    for i, (a, b) in enumerate(zip(spacer.upper(), target20.upper())):
        if a != b:
            score *= POSITION_WEIGHTS[i]
    return score
