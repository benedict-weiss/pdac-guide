"""Characterization tests: the (vectorized) scan must match a naive reference
scanner byte-for-byte, including N bases, lowercase input, and minus-strand
forward-coordinate mapping. Written before vectorization to lock behavior."""

import random

from guide_design.cfd import cfd_score
from guide_design.offtarget import OffTarget, find_offtargets
from guide_design.seq import pam_matches, revcomp

SPACER_LEN, PAM_LEN = 20, 3


def _naive_scan(spacer, frame, strand, pam, max_mismatch):
    out = []
    last = len(frame) - (SPACER_LEN + PAM_LEN)
    for i in range(0, last + 1):
        if not pam_matches(frame[i + SPACER_LEN:i + SPACER_LEN + PAM_LEN], pam):
            continue
        site = frame[i:i + SPACER_LEN]
        mm = sum(1 for x, y in zip(spacer, site) if x != y)
        if mm <= max_mismatch:
            position = i if strand == "+" else last - i
            out.append(OffTarget(strand, position, site, mm, cfd_score(spacer, site)))
    return out


def _naive_find(spacer, reference, pam="NGG", max_mismatch=4):
    ref = reference.upper()
    sp = spacer.upper()
    return (_naive_scan(sp, ref, "+", pam, max_mismatch)
            + _naive_scan(sp, revcomp(ref), "-", pam, max_mismatch))


def _key(hits):
    return [(h.strand, h.position, h.site_seq, h.mismatches, round(h.cfd, 12)) for h in hits]


def test_matches_naive_on_random_references():
    rng = random.Random(1234)
    alphabet = "ACGTN acgt"  # include N, lowercase, and a stray space-free mix
    for _ in range(200):
        n = rng.randint(0, 60)
        reference = "".join(rng.choice("ACGTNacgtn") for _ in range(n))
        spacer = "".join(rng.choice("ACGT") for _ in range(SPACER_LEN))
        mm = rng.randint(0, 6)
        got = find_offtargets(spacer, reference, max_mismatch=mm)
        want = _naive_find(spacer, reference, max_mismatch=mm)
        assert _key(got) == _key(want), f"mismatch on ref={reference!r} mm={mm}"


def test_short_reference_returns_nothing():
    assert find_offtargets("A" * 20, "ACGT") == []
    assert find_offtargets("A" * 20, "") == []


def test_chunk_boundary_matches_naive(monkeypatch):
    import guide_design.offtarget as ot
    monkeypatch.setattr(ot, "_CHUNK", 7)  # force many chunk boundaries
    rng = random.Random(99)
    reference = "".join(rng.choice("ACGTN") for _ in range(300))
    spacer = "".join(rng.choice("ACGT") for _ in range(SPACER_LEN))
    got = ot.find_offtargets(spacer, reference, max_mismatch=5)
    want = _naive_find(spacer, reference, max_mismatch=5)
    assert _key(got) == _key(want)


def test_n_in_reference_counts_as_mismatch():
    spacer = "ACGTACGTACGTACGTACGT"
    # one N inside the protospacer -> exactly 1 mismatch
    site_with_n = "ACGTACGTACGTACGTACGN"
    reference = "AA" + site_with_n + "AGG"
    hits = [h for h in find_offtargets(spacer, reference, max_mismatch=4)
            if h.site_seq == site_with_n]
    assert len(hits) == 1
    assert hits[0].mismatches == 1
