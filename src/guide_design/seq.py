"""Minimal nucleotide-sequence helpers (pure Python, no deps)."""

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    """Reverse complement of a DNA string."""
    return s.translate(_COMPLEMENT)[::-1]


def pam_matches(seq3: str, pattern: str = "NGG") -> bool:
    """True if `seq3` matches `pattern` positionally; `N` is a wildcard."""
    if len(seq3) != len(pattern):
        return False
    return all(p == "N" or p == b for p, b in zip(pattern, seq3.upper()))
