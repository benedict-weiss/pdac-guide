"""Genome / reference off-target search and safety aggregation.

Pure-Python scan suitable for a local locus or a single chromosome (cached).
For genome-scale runs prefer Cas-OFFinder if installed; this module is the
dependency-free fallback. See scripts/fetch_data.py.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cfd import cfd_score
from .seq import revcomp

SPACER_LEN = 20
PAM_LEN = 3
WINDOW = SPACER_LEN + PAM_LEN
# Windows processed per batch; caps the transient mismatch mask at ~CHUNK*20 bytes.
_CHUNK = 2_000_000


@dataclass(frozen=True)
class OffTarget:
    strand: str
    # 0-based start of the matched 23-mer (20-nt protospacer + 3-nt PAM) in the
    # FORWARD reference, reported consistently across both strands.
    position: int
    site_seq: str
    mismatches: int
    cfd: float


def _scan(spacer: str, frame: str, strand: str, pam: str, max_mismatch: int) -> list[OffTarget]:
    last = len(frame) - WINDOW
    # A PAM whose length isn't 3 can never match a 3-nt frame slice (the old
    # pam_matches did a length check), so there are no hits.
    if last < 0 or len(pam) != PAM_LEN:
        return []

    # Encode frame and spacer as ASCII codes. Every base is compared by equality,
    # so N (and any other byte) differs from A/C/G/T and counts as a mismatch,
    # matching the original char-by-char behavior.
    buf = np.frombuffer(frame.encode("ascii"), dtype=np.uint8)
    windows = np.lib.stride_tricks.sliding_window_view(buf, WINDOW)  # (last+1, 23), a view
    spacer_codes = np.frombuffer(spacer.encode("ascii"), dtype=np.uint8)
    pam_codes = pam.upper().encode("ascii")

    # Scan in chunks so the transient boolean mask stays bounded (~CHUNK*20 bytes)
    # even for a whole chromosome, instead of allocating one array over all windows.
    out: list[OffTarget] = []
    n_windows = windows.shape[0]
    for start in range(0, n_windows, _CHUNK):
        block = windows[start : start + _CHUNK]
        mism = (block[:, :SPACER_LEN] != spacer_codes).sum(axis=1)
        pam_ok = np.ones(block.shape[0], dtype=bool)
        for j, code in enumerate(pam_codes):
            if chr(code) != "N":  # 'N' is a wildcard (matches pam_matches)
                pam_ok &= block[:, SPACER_LEN + j] == code
        for local in np.nonzero(pam_ok & (mism <= max_mismatch))[0].tolist():
            i = start + local
            site = frame[i : i + SPACER_LEN]
            # `position` is the 23-mer start in FORWARD reference coordinates. For
            # the minus strand `frame` is revcomp(reference), so map `i` to forward.
            position = i if strand == "+" else last - i
            out.append(OffTarget(strand, position, site, int(mism[local]), cfd_score(spacer, site)))
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
