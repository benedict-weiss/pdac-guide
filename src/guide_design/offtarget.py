"""Genome / reference off-target search and safety aggregation.

Pure-Python scan suitable for a local locus or a single chromosome (cached).
For genome-scale runs prefer Cas-OFFinder if installed; this module is the
dependency-free fallback. See scripts/fetch_data.py.
"""

from __future__ import annotations

from dataclasses import dataclass

from .cfd import cfd_score
from .seq import pam_matches, revcomp

SPACER_LEN = 20
PAM_LEN = 3


@dataclass(frozen=True)
class OffTarget:
    strand: str
    position: int
    site_seq: str
    mismatches: int
    cfd: float


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def _scan(spacer: str, frame: str, strand: str, pam: str, max_mismatch: int) -> list[OffTarget]:
    out: list[OffTarget] = []
    last = len(frame) - (SPACER_LEN + PAM_LEN)
    for i in range(0, last + 1):
        if not pam_matches(frame[i + SPACER_LEN : i + SPACER_LEN + PAM_LEN], pam):
            continue
        site = frame[i : i + SPACER_LEN]
        mm = _hamming(spacer, site)
        if mm <= max_mismatch:
            out.append(OffTarget(strand, i, site, mm, cfd_score(spacer, site)))
    return out


def find_offtargets(
    spacer: str, reference: str, pam: str = "NGG", max_mismatch: int = 4
) -> list[OffTarget]:
    ref = reference.upper()
    plus = _scan(spacer.upper(), ref, "+", pam, max_mismatch)
    minus = _scan(spacer.upper(), revcomp(ref), "-", pam, max_mismatch)
    return plus + minus


def specificity_score(off_cfds: list[float]) -> float:
    """MIT-style aggregate: 100 / (1 + sum of off-target CFDs). 100 = perfectly specific."""
    return 100.0 / (1.0 + sum(off_cfds))


def passes_safety_veto(offtargets: list[OffTarget], cfd_threshold: float = 0.5) -> bool:
    """Fail if any off-target's CFD exceeds the threshold."""
    return all(o.cfd <= cfd_threshold for o in offtargets)
