# PDAC CRISPR Guide-Design Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible Python pipeline that designs and ranks mutation-specific CRISPR guides for the KRAS codon-12 allelic series (G12D/G12V/G12R) in PDAC, optimizing mutant-vs-wild-type discrimination first, then genome-wide off-target safety, then on-target activity.

**Architecture:** A small `guide_design` package decomposed by responsibility: sequence helpers → target definitions → protospacer enumeration → CFD-style scoring → off-target search → ranking/coverage. A thin `scripts/run_analysis.py` wires them into an end-to-end run that emits a recommendations table and figures. Heavy genome-wide off-target search is an optional, cached step (`scripts/fetch_data.py`); the core pipeline runs in minutes on bundled data.

**Tech Stack:** Python ≥3.11, `pandas` (tables), `matplotlib` (figures), `pytest` (tests). Pure-Python sequence logic — no Biopython/numpy required. Optional external tool: Cas-OFFinder (used only if already installed).

## Global Constraints

- **Python ≥ 3.11** (uses `X | None` unions and `@dataclass`).
- **Runtime dependencies:** `pandas>=2.0`, `matplotlib>=3.7` only. Dev: `pytest>=7.0`.
- **No network access at runtime** for the core pipeline (`run_analysis.py`). Network is confined to the optional `scripts/fetch_data.py`.
- **Reproducible one-command run:** `pip install -e ".[dev]" && python scripts/run_analysis.py` must succeed using only bundled data.
- **Cite data sources** in `data/*.tsv` headers and the writeup (GRCh38, COSMIC/TCGA-PAAD/GENIE, Doench et al. 2016 for CFD).
- **Honest-limitations framing:** the CFD model is a documented *approximation*; on-target predictors are *cutting*-trained and used only as a tie-break. Never present scores as validated truth.
- **Position convention (used everywhere):** within a 20-nt spacer, index `0` = 5′ end = **PAM-distal**; index `19` = 3′ end = **PAM-proximal** (adjacent to PAM). "Seed" = PAM-proximal half.

---

### Task 1: Project scaffold & packaging

**Files:**
- Create: `pyproject.toml`
- Create: `src/guide_design/__init__.py`
- Create: `tests/test_smoke.py`
- Create: `.gitignore`

**Interfaces:**
- Consumes: nothing.
- Produces: installable package `guide_design` with `__version__: str`.

- [ ] **Step 1: Write `.gitignore`**

```gitignore
__pycache__/
*.pyc
.pytest_cache/
*.egg-info/
build/
dist/
results/
data/offtarget_cache/*.json
.DS_Store
```

- [ ] **Step 2: Write `pyproject.toml`**

```toml
[project]
name = "guide-design"
version = "0.1.0"
description = "Mutation-specific CRISPR guide design for PDAC (KRAS codon 12)"
requires-python = ">=3.11"
dependencies = ["pandas>=2.0", "matplotlib>=3.7"]

[project.optional-dependencies]
dev = ["pytest>=7.0"]

[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
```

- [ ] **Step 3: Write `src/guide_design/__init__.py`**

```python
"""Mutation-specific CRISPR guide design for PDAC (KRAS codon 12)."""

__version__ = "0.1.0"
```

- [ ] **Step 4: Write the smoke test `tests/test_smoke.py`**

```python
import guide_design


def test_package_imports_with_version():
    assert guide_design.__version__ == "0.1.0"
```

- [ ] **Step 5: Install and run the smoke test**

Run: `pip install -e ".[dev]" && pytest tests/test_smoke.py -v`
Expected: PASS (1 passed).

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml src/guide_design/__init__.py tests/test_smoke.py .gitignore
git commit -m "chore: scaffold guide_design package"
```

---

### Task 2: Sequence helpers

**Files:**
- Create: `src/guide_design/seq.py`
- Test: `tests/test_seq.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `revcomp(s: str) -> str`
  - `pam_matches(seq3: str, pattern: str = "NGG") -> bool` — `N` matches any base; comparison is positional and length-checked.

- [ ] **Step 1: Write the failing test `tests/test_seq.py`**

```python
from guide_design.seq import revcomp, pam_matches


def test_revcomp_basic():
    assert revcomp("ATGC") == "GCAT"


def test_revcomp_is_involution():
    s = "ACGTACGTTTGGCCAA"
    assert revcomp(revcomp(s)) == s


def test_pam_matches_ngg():
    assert pam_matches("TGG") is True
    assert pam_matches("AGG") is True
    assert pam_matches("TAG") is False   # 2nd base not G


def test_pam_matches_length_guard():
    assert pam_matches("GG") is False     # wrong length
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_seq.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.seq'`.

- [ ] **Step 3: Write `src/guide_design/seq.py`**

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_seq.py -v`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add src/guide_design/seq.py tests/test_seq.py
git commit -m "feat: add revcomp and PAM-match sequence helpers"
```

---

### Task 3: Target definitions & PDAC KRAS spectrum

**Files:**
- Create: `data/pdac_kras_spectrum.tsv`
- Create: `src/guide_design/targets.py`
- Test: `tests/test_targets.py`

**Interfaces:**
- Consumes: `pdac_kras_spectrum.tsv`.
- Produces:
  - `KRAS_REF_CONTEXT: str` — 90-nt KRAS CDS coding-sense context (CDS nt 1–90), codon 12 at 0-based indices 33–35.
  - `Target` dataclass: `target_id, gene, allele, ref_seq, snv_index, ref_base, alt_base, frequency`; property `mut_seq -> str`.
  - `load_targets(spectrum_path: str | Path = ...) -> list[Target]`
  - `TARGETS: list[Target]` — default-loaded list.

**Background for the implementer:** KRAS codon 12 = `GGT` (Gly). On the coding strand, CDS nt 34/35/36. Mutations are single-base: **G12R** = c.34 G>C (0-based index 33), **G12D** = c.35 G>A (index 34), **G12V** = c.35 G>T (index 34). The reference context below is CDS nt 1–90 of canonical KRAS; verify `KRAS_REF_CONTEXT[33:36] == "GGT"`.

- [ ] **Step 1: Write `data/pdac_kras_spectrum.tsv`**

```tsv
# KRAS codon-12 mutation spectrum in PDAC. freq_of_kras_mut = approx fraction
# of KRAS-mutant PDAC tumors carrying this allele. Sources: COSMIC, TCGA-PAAD,
# AACR GENIE (approximate; see writeup for citations). snv_index is 0-based into
# KRAS_REF_CONTEXT (CDS nt 1-90, coding sense).
allele	cds_pos	ref_base	alt_base	snv_index	freq_of_kras_mut
G12D	35	G	A	34	0.40
G12V	35	G	T	34	0.30
G12R	34	G	C	33	0.17
```

- [ ] **Step 2: Write the failing test `tests/test_targets.py`**

```python
from guide_design.targets import KRAS_REF_CONTEXT, Target, load_targets, TARGETS


def test_reference_context_codon12_is_ggt():
    # Codon 12 sits at 0-based indices 33-35 of the CDS context.
    assert len(KRAS_REF_CONTEXT) == 90
    assert KRAS_REF_CONTEXT[33:36] == "GGT"


def test_targets_loaded_for_three_alleles():
    alleles = {t.allele for t in TARGETS}
    assert alleles == {"G12D", "G12V", "G12R"}


def test_ref_base_matches_context():
    for t in TARGETS:
        assert KRAS_REF_CONTEXT[t.snv_index] == t.ref_base


def test_mut_seq_applies_single_substitution():
    g12d = next(t for t in TARGETS if t.allele == "G12D")
    mut = g12d.mut_seq
    assert mut[34] == "A"                 # G>A applied
    assert mut[:34] == KRAS_REF_CONTEXT[:34]   # nothing else changed
    assert mut[35:] == KRAS_REF_CONTEXT[35:]
    assert mut[33:36] == "GAT"            # codon 12 now Asp


def test_g12r_changes_first_codon_base():
    g12r = next(t for t in TARGETS if t.allele == "G12R")
    assert g12r.snv_index == 33
    assert g12r.mut_seq[33:36] == "CGT"   # Arg
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pytest tests/test_targets.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.targets'`.

- [ ] **Step 4: Write `src/guide_design/targets.py`**

```python
"""KRAS codon-12 targets and the PDAC mutation spectrum."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# Canonical KRAS CDS, nucleotides 1-90 (coding sense, 5'->3'). Codon 12 (GGT)
# is at 0-based indices 33-35. Source: GRCh38 / canonical KRAS transcript.
KRAS_REF_CONTEXT: str = (
    "ATGACTGAAT"  # 1-10
    "ATAAACTTGT"  # 11-20
    "GGTAGTTGGA"  # 21-30
    "GCTGGTGGCG"  # 31-40  (codon 12 GGT at idx 33-35)
    "TAGGCAAGAG"  # 41-50
    "TGCCTTGACG"  # 51-60
    "ATACAGCTAA"  # 61-70
    "TTCAGAATCA"  # 71-80
    "TTTTGTGGAC"  # 81-90
)

_DEFAULT_SPECTRUM = Path(__file__).resolve().parents[2] / "data" / "pdac_kras_spectrum.tsv"


@dataclass(frozen=True)
class Target:
    target_id: str
    gene: str
    allele: str
    ref_seq: str
    snv_index: int       # 0-based index into ref_seq of the changed base
    ref_base: str
    alt_base: str
    frequency: float     # fraction of KRAS-mutant PDAC tumors

    @property
    def mut_seq(self) -> str:
        """Reference context with the single allele substitution applied."""
        i = self.snv_index
        return self.ref_seq[:i] + self.alt_base + self.ref_seq[i + 1:]


def load_targets(spectrum_path: str | Path = _DEFAULT_SPECTRUM) -> list[Target]:
    df = pd.read_csv(spectrum_path, sep="\t", comment="#")
    targets: list[Target] = []
    for row in df.itertuples(index=False):
        targets.append(
            Target(
                target_id=f"KRAS_{row.allele}",
                gene="KRAS",
                allele=row.allele,
                ref_seq=KRAS_REF_CONTEXT,
                snv_index=int(row.snv_index),
                ref_base=str(row.ref_base),
                alt_base=str(row.alt_base),
                frequency=float(row.freq_of_kras_mut),
            )
        )
    return targets


TARGETS: list[Target] = load_targets()
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pytest tests/test_targets.py -v`
Expected: PASS (5 passed).

- [ ] **Step 6: Commit**

```bash
git add data/pdac_kras_spectrum.tsv src/guide_design/targets.py tests/test_targets.py
git commit -m "feat: define KRAS codon-12 targets and PDAC spectrum"
```

---

### Task 4: Protospacer enumeration

**Files:**
- Create: `src/guide_design/enumerate.py`
- Test: `tests/test_enumerate.py`

**Interfaces:**
- Consumes: `Target` (Task 3), `revcomp`/`pam_matches` (Task 2).
- Produces:
  - `Guide` dataclass: `target_id, allele, strand ('+'|'-'), start (int, index in working frame), spacer (20-nt, == mutant protospacer), pam (3-nt mutant PAM), wt_protospacer (20-nt WT at same site, spacer orientation), wt_has_pam (bool), snv_pos (int 1..20 PAM-distal→proximal, or None if SNV is in the PAM only), mutant_specific_pam (bool)`.
  - `enumerate_guides(target: Target, pam: str = "NGG") -> list[Guide]`

**Design notes for the implementer:**
- Enumerate over the **mutant** sequence (so the PAM reflects the mutant). A candidate is kept only if the SNV lies in the protospacer **or** in the PAM — otherwise it cannot discriminate mutant from WT.
- `spacer` equals the mutant protospacer (the gRNA is complementary to it by construction), so on-mutant recognition is a perfect match.
- For the minus strand, reverse-complement both the mutant and WT frames and re-run the same forward scan; map the SNV index into the reversed frame as `len-1-snv_index`.
- `mutant_specific_pam = not wt_has_pam` (mutant always has a PAM by construction, since we enumerate on the mutant). This is the gold discrimination class.

- [ ] **Step 1: Write the failing test `tests/test_enumerate.py`**

```python
from guide_design.targets import Target
from guide_design.enumerate import Guide, enumerate_guides


def _target(ref_seq, snv_index, ref_base, alt_base, allele="X"):
    return Target("T_" + allele, "GENE", allele, ref_seq, snv_index, ref_base, alt_base, 1.0)


def test_plus_strand_guide_with_snv_in_seed():
    # 20 A's, then PAM 'TGG'. SNV at index 19 (PAM-proximal end of protospacer).
    ref = "A" * 20 + "TGGAAA"
    t = _target(ref, 19, "A", "G")
    guides = enumerate_guides(t)
    g = next(x for x in guides if x.strand == "+" and x.start == 0)
    assert g.pam == "TGG"
    assert g.spacer == "A" * 19 + "G"       # mutant protospacer
    assert g.wt_protospacer == "A" * 20
    assert g.wt_has_pam is True
    assert g.snv_pos == 20                   # PAM-proximal (seed)
    assert g.mutant_specific_pam is False


def test_mutant_specific_pam_when_snv_creates_pam():
    # PAM region is 'TAG' in WT (no NGG); SNV A>G at index 21 makes 'TGG'.
    ref = "A" * 20 + "TAGAAA"
    t = _target(ref, 21, "A", "G")
    guides = enumerate_guides(t)
    g = next(x for x in guides if x.strand == "+" and x.start == 0)
    assert g.pam == "TGG"
    assert g.spacer == "A" * 20
    assert g.wt_has_pam is False             # WT 'TAG' is not a PAM
    assert g.mutant_specific_pam is True
    assert g.snv_pos is None                 # SNV is in the PAM, not the protospacer


def test_guides_not_covering_snv_are_excluded():
    # SNV far from any protospacer/PAM overlap -> no guide should reference it.
    ref = "TGG" + "A" * 30      # only PAM is at the very start, no 20-mer before it
    t = _target(ref, 25, "A", "G")
    guides = enumerate_guides(t)
    assert guides == []


def test_minus_strand_guide_found():
    # 'CCA' on the forward strand => 'TGG' PAM on the reverse strand.
    ref = "CCA" + "T" * 20 + "GGG"
    t = _target(ref, 4, "T", "C")   # SNV inside the reverse-strand protospacer
    guides = enumerate_guides(t)
    assert any(g.strand == "-" for g in guides)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_enumerate.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.enumerate'`.

- [ ] **Step 3: Write `src/guide_design/enumerate.py`**

```python
"""Enumerate candidate protospacers that discriminate a mutation from WT."""

from __future__ import annotations

from dataclasses import dataclass

from .seq import pam_matches, revcomp
from .targets import Target

SPACER_LEN = 20
PAM_LEN = 3


@dataclass(frozen=True)
class Guide:
    target_id: str
    allele: str
    strand: str            # '+' or '-'
    start: int             # protospacer start index within the working frame
    spacer: str            # 20-nt, equals the mutant protospacer (gRNA target)
    pam: str               # 3-nt mutant PAM
    wt_protospacer: str    # 20-nt WT sequence at the same site, spacer orientation
    wt_has_pam: bool
    snv_pos: int | None    # 1..20 (PAM-distal->proximal) if SNV in protospacer, else None
    mutant_specific_pam: bool


def _enumerate_strand(
    target: Target, mut_frame: str, wt_frame: str, snv_index: int, strand: str, pam: str
) -> list[Guide]:
    guides: list[Guide] = []
    last = len(mut_frame) - (SPACER_LEN + PAM_LEN)
    for i in range(0, last + 1):
        mut_pam = mut_frame[i + SPACER_LEN : i + SPACER_LEN + PAM_LEN]
        if not pam_matches(mut_pam, pam):
            continue
        snv_in_proto = i <= snv_index <= i + SPACER_LEN - 1
        snv_in_pam = i + SPACER_LEN <= snv_index <= i + SPACER_LEN + PAM_LEN - 1
        if not (snv_in_proto or snv_in_pam):
            continue
        wt_pam = wt_frame[i + SPACER_LEN : i + SPACER_LEN + PAM_LEN]
        wt_has_pam = pam_matches(wt_pam, pam)
        snv_pos = (snv_index - i) + 1 if snv_in_proto else None
        guides.append(
            Guide(
                target_id=target.target_id,
                allele=target.allele,
                strand=strand,
                start=i,
                spacer=mut_frame[i : i + SPACER_LEN],
                pam=mut_pam,
                wt_protospacer=wt_frame[i : i + SPACER_LEN],
                wt_has_pam=wt_has_pam,
                snv_pos=snv_pos,
                mutant_specific_pam=not wt_has_pam,
            )
        )
    return guides


def enumerate_guides(target: Target, pam: str = "NGG") -> list[Guide]:
    mut, wt = target.mut_seq, target.ref_seq
    plus = _enumerate_strand(target, mut, wt, target.snv_index, "+", pam)
    rc_index = len(mut) - 1 - target.snv_index
    minus = _enumerate_strand(target, revcomp(mut), revcomp(wt), rc_index, "-", pam)
    return plus + minus
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_enumerate.py -v`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add src/guide_design/enumerate.py tests/test_enumerate.py
git commit -m "feat: enumerate discriminating protospacers on both strands"
```

---

### Task 5: CFD-style scoring & discrimination

**Files:**
- Create: `src/guide_design/cfd.py`
- Create: `src/guide_design/score.py`
- Test: `tests/test_score.py`

**Interfaces:**
- Consumes: `Guide` (Task 4).
- Produces:
  - `cfd.POSITION_WEIGHTS: list[float]` (len 20; index 0 = PAM-distal), `cfd.cfd_score(spacer: str, target20: str) -> float`.
  - `score.DiscriminationResult` dataclass: `score_mut, score_wt, margin, mechanism`.
  - `score.discrimination(guide: Guide) -> DiscriminationResult`
  - `score.on_target_proxy(spacer: str) -> float`
  - `score.SEED_START: int = 11`

**Design notes:** `cfd_score` multiplies a per-position "activity retained" weight for every mismatched position (perfect match → 1.0). The weight vector is a documented monotonic approximation of the Doench et al. 2016 CFD averages (PAM-distal tolerant → PAM-proximal/seed intolerant). It is explicitly an approximation; the writeup must say so. `discrimination` sets `score_mut = 1.0` (gRNA matches the mutant by construction), and `score_wt = 0.0` when the WT lacks a PAM (mutant-specific PAM — WT cannot be bound), else `cfd_score(spacer, wt_protospacer)`.

- [ ] **Step 1: Write the failing test `tests/test_score.py`**

```python
from guide_design.targets import Target
from guide_design.enumerate import enumerate_guides
from guide_design.cfd import cfd_score, POSITION_WEIGHTS
from guide_design.score import discrimination, on_target_proxy


def test_cfd_perfect_match_is_one():
    assert cfd_score("A" * 20, "A" * 20) == 1.0


def test_cfd_seed_mismatch_lower_than_distal():
    spacer = "A" * 20
    distal = "C" + "A" * 19    # mismatch at index 0 (PAM-distal)
    seed = "A" * 19 + "C"      # mismatch at index 19 (PAM-proximal)
    assert cfd_score(spacer, distal) > cfd_score(spacer, seed)


def test_position_weights_monotone_nonincreasing():
    assert len(POSITION_WEIGHTS) == 20
    assert all(POSITION_WEIGHTS[i] >= POSITION_WEIGHTS[i + 1] for i in range(19))


def _target(ref, idx, rb, ab, allele="X"):
    return Target("T", "G", allele, ref, idx, rb, ab, 1.0)


def test_discrimination_seed_mechanism_high_margin():
    # SNV at PAM-proximal end -> seed mismatch -> high margin.
    t = _target("A" * 20 + "TGGAAA", 19, "A", "G")
    g = next(x for x in enumerate_guides(t) if x.strand == "+" and x.snv_pos == 20)
    d = discrimination(g)
    assert d.score_mut == 1.0
    assert d.margin > 0.9
    assert d.mechanism == "seed_mismatch"


def test_discrimination_mutant_specific_pam_zero_wt():
    t = _target("A" * 20 + "TAGAAA", 21, "A", "G")
    g = next(x for x in enumerate_guides(t) if x.mutant_specific_pam)
    d = discrimination(g)
    assert d.score_wt == 0.0
    assert d.margin == 1.0
    assert d.mechanism == "pam_differential"


def test_on_target_proxy_penalizes_polyt():
    assert on_target_proxy("AAGCTAGCTAGCTAGCTAGC") > on_target_proxy("TTTTAGCTAGCTAGCTAGCT")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_score.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.cfd'`.

- [ ] **Step 3: Write `src/guide_design/cfd.py`**

```python
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
```

- [ ] **Step 4: Write `src/guide_design/score.py`**

```python
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
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_score.py -v`
Expected: PASS (6 passed).

- [ ] **Step 6: Commit**

```bash
git add src/guide_design/cfd.py src/guide_design/score.py tests/test_score.py
git commit -m "feat: add CFD-style scoring and mutant-vs-WT discrimination"
```

---

### Task 6: Off-target search & safety veto

**Files:**
- Create: `src/guide_design/offtarget.py`
- Test: `tests/test_offtarget.py`

**Interfaces:**
- Consumes: `cfd_score` (Task 5), `pam_matches`/`revcomp` (Task 2).
- Produces:
  - `OffTarget` dataclass: `strand, position, site_seq, mismatches, cfd`.
  - `find_offtargets(spacer: str, reference: str, pam: str = "NGG", max_mismatch: int = 4) -> list[OffTarget]`
  - `specificity_score(off_cfds: list[float]) -> float` — MIT-style aggregate in (0, 100].
  - `passes_safety_veto(offtargets: list[OffTarget], cfd_threshold: float = 0.5) -> bool`

**Design notes:** scan both strands of `reference` for 20-nt sites that carry a valid PAM and lie within `max_mismatch` of the spacer; score each with `cfd_score`. `position`/`strand` locate the site; `mismatches` is the Hamming distance. The on-target site appears as a 0-mismatch hit — callers exclude `mismatches == 0`. The veto fails if any off-target's CFD exceeds `cfd_threshold`. The **WT KRAS allele** is evaluated separately by the discrimination step (Task 5) and surfaced in reporting (Task 8); the genome scan here catches sites elsewhere.

- [ ] **Step 1: Write the failing test `tests/test_offtarget.py`**

```python
from guide_design.offtarget import (
    OffTarget, find_offtargets, specificity_score, passes_safety_veto,
)


def test_finds_exact_on_target_site():
    spacer = "ACGTACGTACGTACGTACGT"
    reference = "TTT" + spacer + "TGG" + "TTT"   # spacer + NGG PAM
    hits = find_offtargets(spacer, reference, max_mismatch=4)
    exact = [h for h in hits if h.mismatches == 0]
    assert len(exact) == 1
    assert exact[0].strand == "+"


def test_finds_near_match_within_threshold_only():
    spacer = "ACGTACGTACGTACGTACGT"
    near = "ACGTACGTACGTACGTACGA"          # 1 mismatch at PAM-proximal end
    far = "TTTTTTTTTTTTTTTTTTTT"           # >4 mismatches
    reference = "AA" + near + "AGG" + "AA" + far + "CGG"
    hits = find_offtargets(spacer, reference, max_mismatch=4)
    seqs = {h.site_seq for h in hits}
    assert near in seqs
    assert far not in seqs


def test_specificity_score_decreases_with_offtargets():
    assert specificity_score([]) == 100.0
    assert specificity_score([0.1]) < 100.0
    assert specificity_score([0.8]) < specificity_score([0.1])


def test_safety_veto_fails_on_high_cfd_offtarget():
    good = [OffTarget("+", 10, "A" * 20, 3, 0.2)]
    bad = [OffTarget("+", 10, "A" * 20, 1, 0.8)]
    assert passes_safety_veto(good, cfd_threshold=0.5) is True
    assert passes_safety_veto(bad, cfd_threshold=0.5) is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_offtarget.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.offtarget'`.

- [ ] **Step 3: Write `src/guide_design/offtarget.py`**

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_offtarget.py -v`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
git add src/guide_design/offtarget.py tests/test_offtarget.py
git commit -m "feat: add off-target search, specificity score, safety veto"
```

---

### Task 7: Ranking & patient coverage

**Files:**
- Create: `src/guide_design/rank.py`
- Create: `src/guide_design/coverage.py`
- Test: `tests/test_rank.py`
- Test: `tests/test_coverage.py`

**Interfaces:**
- Consumes: `Guide` (Task 4), `discrimination`/`on_target_proxy` (Task 5), `find_offtargets`/`specificity_score`/`passes_safety_veto` (Task 6), `Target` (Task 3).
- Produces:
  - `rank.RankedGuide` dataclass: `guide, margin, score_wt, mechanism, specificity (float|None), on_target, veto_pass`.
  - `rank.evaluate_guide(guide, reference: str | None = None, max_mismatch=4, cfd_threshold=0.5) -> RankedGuide`
  - `rank.rank_guides(rgs: list[RankedGuide], margin_threshold: float = 0.5) -> list[RankedGuide]`
  - `rank.recommend_per_allele(rgs: list[RankedGuide], margin_threshold: float = 0.5) -> dict[str, RankedGuide | None]`
  - `coverage.coverage_curve(recommendations, targets, kras_fraction: float = 0.90) -> list[dict]`

**Design notes:** ranking is lexicographic and transparent — a guide *qualifies* iff `veto_pass and margin >= margin_threshold`; qualified guides sort ahead of unqualified, then by `(margin, on_target)` descending. `recommend_per_allele` returns the best qualified guide per allele, or `None` if the allele is **undesignable** (honest negative result). When `reference is None`, off-target is not evaluated: `specificity = None`, `veto_pass = True` (so the local demo still ranks by discrimination).

- [ ] **Step 1: Write the failing test `tests/test_rank.py`**

```python
from guide_design.targets import Target
from guide_design.enumerate import enumerate_guides
from guide_design.rank import evaluate_guide, rank_guides, recommend_per_allele


def _target(ref, idx, rb, ab, allele="X", freq=1.0):
    return Target("KRAS_" + allele, "KRAS", allele, ref, idx, rb, ab, freq)


def test_evaluate_guide_without_reference_skips_offtarget():
    t = _target("A" * 20 + "TGGAAA", 19, "A", "G", allele="G12D")
    g = enumerate_guides(t)[0]
    rg = evaluate_guide(g, reference=None)
    assert rg.specificity is None
    assert rg.veto_pass is True
    assert rg.margin > 0.9


def test_rank_orders_qualified_first_then_by_margin():
    t_seed = _target("A" * 20 + "TGGAAA", 19, "A", "G", allele="A")   # seed -> high margin
    t_distal = _target("C" + "A" * 19 + "TGGAAA", 0, "C", "G", allele="B")  # distal -> low margin
    rgs = [evaluate_guide(enumerate_guides(t_distal)[0]),
           evaluate_guide(enumerate_guides(t_seed)[0])]
    ranked = rank_guides(rgs, margin_threshold=0.5)
    assert ranked[0].margin > ranked[1].margin    # high-margin seed guide first


def test_recommend_marks_undesignable_allele_none():
    # Only a distal mismatch available -> margin below threshold -> undesignable.
    t = _target("C" + "A" * 19 + "TGGAAA", 0, "C", "G", allele="G12X")
    rgs = [evaluate_guide(g) for g in enumerate_guides(t)]
    rec = recommend_per_allele(rgs, margin_threshold=0.9)
    assert rec["G12X"] is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_rank.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.rank'`.

- [ ] **Step 3: Write `src/guide_design/rank.py`**

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_rank.py -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Write the failing test `tests/test_coverage.py`**

```python
from guide_design.targets import Target
from guide_design.enumerate import enumerate_guides
from guide_design.rank import evaluate_guide, recommend_per_allele
from guide_design.coverage import coverage_curve


def _designable_target(allele, freq):
    # Seed-position SNV -> a qualifying guide always exists.
    return Target("KRAS_" + allele, "KRAS", allele, "A" * 20 + "TGGAAA", 19, "A", "G", freq)


def test_coverage_accumulates_in_frequency_order():
    targets = [_designable_target("G12D", 0.40),
               _designable_target("G12V", 0.30),
               _designable_target("G12R", 0.17)]
    rgs = [evaluate_guide(g) for t in targets for g in enumerate_guides(t)]
    rec = recommend_per_allele(rgs, margin_threshold=0.5)
    rows = coverage_curve(rec, targets, kras_fraction=0.90)
    assert [r["added_allele"] for r in rows] == ["G12D", "G12V", "G12R"]
    assert abs(rows[-1]["cumulative_kras_frac"] - 0.87) < 1e-9
    assert abs(rows[-1]["cumulative_pdac_frac"] - 0.87 * 0.90) < 1e-9
    assert abs(rows[0]["marginal_pdac"] - 0.40 * 0.90) < 1e-9


def test_coverage_skips_undesignable_alleles():
    targets = [_designable_target("G12D", 0.40)]
    rec = {"G12D": None}
    assert coverage_curve(rec, targets) == []
```

- [ ] **Step 6: Run test to verify it fails**

Run: `pytest tests/test_coverage.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.coverage'`.

- [ ] **Step 7: Write `src/guide_design/coverage.py`**

```python
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
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `pytest tests/test_coverage.py -v`
Expected: PASS (2 passed).

- [ ] **Step 9: Commit**

```bash
git add src/guide_design/rank.py src/guide_design/coverage.py tests/test_rank.py tests/test_coverage.py
git commit -m "feat: add transparent ranking and patient-coverage analysis"
```

---

### Task 8: End-to-end pipeline, figures & optional data fetch

**Files:**
- Create: `src/guide_design/pipeline.py`
- Create: `scripts/run_analysis.py`
- Create: `scripts/fetch_data.py`
- Test: `tests/test_pipeline.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `pipeline.run(targets, reference: str | None, margin_threshold=0.5) -> tuple[dict[str, RankedGuide|None], pd.DataFrame, list[dict]]` returning `(recommendations, table, coverage_rows)`.
  - `scripts/run_analysis.py` CLI → writes `results/recommendations.tsv` + `results/figures/*.png`.
  - `scripts/fetch_data.py` CLI → optional genome/chr12 download + cached off-target search.

**Design notes:** `pipeline.run` enumerates guides for every target, evaluates each, builds a tidy DataFrame (one row per allele with the recommended guide + key scores, plus an `undesignable` flag), and computes the coverage curve. `run_analysis.py` loads bundled `TARGETS`, optionally loads an off-target reference (a FASTA path or the cached JSON; if absent it warns and proceeds with `reference=None`), then writes outputs. `fetch_data.py` documents and performs the heavy step; it is **not** exercised by unit tests (network-bound).

- [ ] **Step 1: Write the failing test `tests/test_pipeline.py`**

```python
import pandas as pd

from guide_design.targets import Target
from guide_design.pipeline import run


def _designable(allele, freq):
    return Target("KRAS_" + allele, "KRAS", allele, "A" * 20 + "TGGAAA", 19, "A", "G", freq)


def test_pipeline_returns_row_per_allele_and_recommends_designable():
    targets = [_designable("G12D", 0.40), _designable("G12V", 0.30)]
    # Small synthetic reference so the off-target path runs deterministically.
    reference = "GGGTTT" + "A" * 20 + "TGG" + "CCCCCC"
    rec, table, cov = run(targets, reference=reference, margin_threshold=0.5)
    assert isinstance(table, pd.DataFrame)
    assert set(table["allele"]) == {"G12D", "G12V"}
    assert (~table["undesignable"]).all()
    assert table["recommended_spacer"].notna().all()
    assert len(cov) == 2


def test_pipeline_flags_undesignable_allele():
    # Distal-only mismatch + high threshold -> undesignable.
    t = Target("KRAS_G12X", "KRAS", "G12X", "C" + "A" * 19 + "TGGAAA", 0, "C", "G", 0.1)
    rec, table, cov = run([t], reference=None, margin_threshold=0.9)
    assert bool(table.loc[table["allele"] == "G12X", "undesignable"].iloc[0]) is True
    assert cov == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_pipeline.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'guide_design.pipeline'`.

- [ ] **Step 3: Write `src/guide_design/pipeline.py`**

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_pipeline.py -v`
Expected: PASS (2 passed).

- [ ] **Step 5: Write `scripts/run_analysis.py`**

```python
"""Run the PDAC guide-design pipeline and write results + figures.

Usage:
    python scripts/run_analysis.py [--reference path/to/offtarget.fa]

Without --reference (and no cached off-target data) the genome-wide off-target
step is skipped and clearly reported; discrimination + coverage still run.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from guide_design.pipeline import run
from guide_design.targets import TARGETS

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"


def _load_reference(path: str | None) -> str | None:
    if not path:
        print("WARNING: no --reference given; genome-wide off-target search SKIPPED. "
              "Run scripts/fetch_data.py to build a cached reference.")
        return None
    seq = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", default=None, help="FASTA for off-target search")
    ap.add_argument("--margin-threshold", type=float, default=0.5)
    args = ap.parse_args()

    reference = _load_reference(args.reference)
    _, table, cov = run(TARGETS, reference=reference, margin_threshold=args.margin_threshold)

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    table.to_csv(RESULTS / "recommendations.tsv", sep="\t", index=False)
    print(table.to_string(index=False))

    # Figure 1: discrimination margin per allele.
    des = table[~table["undesignable"]]
    plt.figure()
    plt.bar(des["allele"], des["margin"])
    plt.ylabel("Discrimination margin (mut - WT)")
    plt.title("Per-allele discrimination margin")
    plt.ylim(0, 1)
    plt.savefig(FIGURES / "discrimination_margin.png", dpi=120, bbox_inches="tight")
    plt.close()

    # Figure 2: cumulative patient coverage.
    if cov:
        labels = [r["added_allele"] for r in cov]
        cum = [r["cumulative_pdac_frac"] for r in cov]
        plt.figure()
        plt.plot(range(1, len(labels) + 1), cum, marker="o")
        plt.xticks(range(1, len(labels) + 1), ["+" + l for l in labels])
        plt.ylabel("Cumulative fraction of PDAC patients")
        plt.title("Patient coverage vs guide set")
        plt.ylim(0, 1)
        plt.savefig(FIGURES / "coverage_curve.png", dpi=120, bbox_inches="tight")
        plt.close()

    if cov:
        print(f"\nCoverage with {len(cov)} guides: "
              f"{cov[-1]['cumulative_pdac_frac']*100:.1f}% of PDAC patients.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 6: Run the pipeline end-to-end**

Run: `python scripts/run_analysis.py`
Expected: prints a table with one row per allele, writes `results/recommendations.tsv` and `results/figures/*.png`, and prints a coverage line. (Off-target columns are blank with the SKIPPED warning — expected without a reference.)

- [ ] **Step 7: Write `scripts/fetch_data.py`**

```python
"""OPTIONAL: build a reference for genome-wide off-target search.

Downloads a reference FASTA (default: human chromosome 12 from Ensembl, which
contains KRAS) so off-target search runs against real sequence. Full-genome is
the documented extension (concatenate all chromosomes, or use Cas-OFFinder).

Usage:
    python scripts/fetch_data.py --out data/chr12.fa
    python scripts/run_analysis.py --reference data/chr12.fa

NOTE: chromosome 12 alone is a deliberate, documented scope limit -- true
off-targets can occur on any chromosome. State this in the writeup.
"""

from __future__ import annotations

import argparse
import urllib.request

# Ensembl REST: soft-masked chromosome sequence. Large download (~tens of MB).
ENSEMBL_CHR12 = (
    "https://rest.ensembl.org/sequence/region/human/12"
    "?content-type=text/x-fasta"
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/chr12.fa")
    ap.add_argument("--url", default=ENSEMBL_CHR12)
    args = ap.parse_args()
    print(f"Downloading {args.url} -> {args.out} (this is the heavy, optional step)...")
    urllib.request.urlretrieve(args.url, args.out)
    print("Done. Re-run: python scripts/run_analysis.py --reference", args.out)


if __name__ == "__main__":
    main()
```

- [ ] **Step 8: Commit**

```bash
git add src/guide_design/pipeline.py scripts/run_analysis.py scripts/fetch_data.py tests/test_pipeline.py
git commit -m "feat: end-to-end pipeline, figures, optional data fetch"
```

---

### Task 9: README & writeup

**Files:**
- Create: `README.md`
- Create: `writeup/writeup.md`

**Interfaces:**
- Consumes: the committed code + a real run's `results/recommendations.tsv`.
- Produces: documentation deliverables. No tests.

**Note for the implementer:** Run `python scripts/run_analysis.py` first and read `results/recommendations.tsv`. Fill every `<...>` placeholder below with the actual computed values (recommended spacers, margins, mechanisms, coverage %). Do not invent numbers — copy them from the results. The reasoning/decision-tree prose is pre-written; only the computed figures need filling.

- [ ] **Step 1: Write `README.md`**

````markdown
# PDAC CRISPR Guide Design (KRAS codon 12)

Mutation-specific guide-RNA design for a conditional CRISPR therapy: recognition
of a cancer mutation triggers a kill switch, so the guide carries both efficacy
and safety. This repo designs and ranks guides for the KRAS codon-12 alleles
(G12D / G12V / G12R) that dominate PDAC.

## Quickstart

```bash
pip install -e ".[dev]"
pytest                       # run the test suite
python scripts/run_analysis.py
```

Outputs: `results/recommendations.tsv` and `results/figures/*.png`.

## Genome-wide off-target search (optional, heavier)

The core run scores discrimination and the wild-type-allele safety check without
any download. For a genome-scale off-target scan:

```bash
python scripts/fetch_data.py --out data/chr12.fa     # downloads chr12 (contains KRAS)
python scripts/run_analysis.py --reference data/chr12.fa
```

Chromosome 12 is a documented scope limit; full-genome (all chromosomes, or
Cas-OFFinder if installed) is the extension.

## How it works (decision tree)

1. **Target = KRAS.** Truncal in ~90% of PDAC, a single recurrent hotspot codon,
   and an essential driver (hard to lose as an escape route).
2. **Per-allele guides.** Each allele is a different DNA sequence; we design one
   guide per allele and report patient coverage.
3. **Discrimination-first scoring.** Enumerate protospacers (both strands, a PAM
   panel); score the mutant-vs-WT margin first; apply a genome-wide off-target
   veto; break ties with an on-target proxy.
4. **Honest judgement.** Report achievable margins, undesignable alleles, the
   WT-allele safety check, and a coverage curve.

## Module map

| File | Responsibility |
|------|----------------|
| `seq.py` | reverse-complement, PAM matching |
| `targets.py` | KRAS targets + PDAC mutation spectrum |
| `enumerate.py` | protospacer enumeration (both strands, SNV positioning) |
| `cfd.py` / `score.py` | CFD-style discrimination scoring + on-target proxy |
| `offtarget.py` | genome/reference off-target search + safety veto |
| `rank.py` / `coverage.py` | ranking + patient coverage |
| `pipeline.py` + `scripts/` | end-to-end run, figures, optional fetch |

## Limitations

See `writeup/writeup.md`. In short: the CFD model is a documented approximation;
on-target predictors are cutting-trained while the trigger is binding; locus
accessibility, WT-KRAS LOH, delivery, and wet-lab validation are out of scope.
````

- [ ] **Step 2: Write `writeup/writeup.md`** (≤2 pages; fill bracketed values from `results/recommendations.tsv`)

```markdown
# Designing a mutation-specific CRISPR therapy for PDAC

## The problem, reframed
FinalDose's molecule fires its payload *on recognition*, so the guide RNA is the
whole therapy: efficacy and safety both live in it. That turns "guide design"
into three coupled questions — what to target, how to tell the mutant from the
wild-type (WT) allele one base away, and how to avoid firing anywhere else in the
genome. The deadliest off-target is not random: it is the WT allele present in
every healthy cell.

## Target selection: why KRAS codon 12
PDAC is ~90% KRAS-mutant, concentrated at a single hotspot codon (12). KRAS is an
essential driver, so the tumor cannot easily delete it to escape — unlike TP53,
CDKN2A, or SMAD4, which are frequently lost by deletion or are subclonal and make
poor point-mutation triggers. So we design at codon 12 and argue the rest in
prose. Within codon 12, the common alleles are G12D (~40%), G12V (~30%), and
G12R (~17%) of KRAS-mutant tumors; each is a distinct DNA sequence needing its
own guide.

## Method: discrimination first
For each allele we enumerate every candidate protospacer (both strands, NGG PAM,
extensible to other PAMs) whose recognition overlaps the mutation, then score:
1. **Discrimination margin** Δ = score(mutant) − score(WT), using a CFD-style
   position-weighted mismatch model. A mismatch in the PAM-proximal seed, or a
   mutation that gives the mutant a PAM the WT lacks, drives WT recognition to ~0.
2. **Off-target safety veto** — a genome/chr12 scan, CFD-aggregated, vetoing any
   guide with a high-CFD off-target; the WT allele is checked explicitly.
3. **On-target proxy** — a documented guardrail (GC, homopolymers, Pol III
   terminator), used only to break ties. It is *not* a validated predictor.

Ranking is lexicographic and inspectable, not a black box.

## Recommendation
| Allele | Spacer (5'→3') | Strand | Mechanism | Margin | Coverage added |
|--------|----------------|--------|-----------|--------|----------------|
| G12D | `<from results>` | `<+/->` | `<mechanism>` | `<margin>` | `<marginal_pdac>` |
| G12V | `<from results>` | `<+/->` | `<mechanism>` | `<margin>` | `<marginal_pdac>` |
| G12R | `<from results>` | `<+/->` | `<mechanism>` | `<margin>` | `<marginal_pdac>` |

The three-guide set addresses ~`<cumulative_pdac_frac>`% of PDAC patients.
`<Note any allele flagged undesignable and why — e.g. no PAM placing the SNV in
the seed.>`

## Trade-offs I weighed
- **Discrimination vs coverage:** chasing rarer alleles adds patients but often
  worse margins; I prioritized clean discrimination over coverage.
- **Seed mismatch vs mutant-specific PAM:** PAM-differential guides give the
  cleanest separation but are not always available per allele.
- **Reproducibility vs scale:** the default run needs no download; genome-scale
  off-target is an explicit, cached extension.

## Honest limitations / with more time
- **Model:** CFD is an approximation of *cutting* data; the trigger is *binding*.
  I would substitute the full empirical CFD/CRISPRoff-style tables and a
  binding-specific model, and validate in isogenic mutant-vs-WT reporter lines.
- **WT-KRAS LOH / zygosity:** if a tumor loses the WT allele, discrimination
  matters less; copy-number context should feed target choice.
- **Locus accessibility & expression**, guide secondary structure, delivery/PK,
  and immunogenicity are out of scope here.
- **Off-target scope:** chr12 by default; true safety needs all chromosomes plus
  common population variants (a guide can be specific in the reference but not in
  a patient whose WT carries a SNP).
```

- [ ] **Step 3: Run the full suite and the pipeline once more to confirm green**

Run: `pytest -q && python scripts/run_analysis.py`
Expected: all tests pass; results + figures regenerate; coverage line prints.

- [ ] **Step 4: Commit**

```bash
git add README.md writeup/writeup.md results/recommendations.tsv
git commit -m "docs: add README and writeup with computed recommendations"
```

---

## Self-Review

**1. Spec coverage**

| Spec element | Task |
|--------------|------|
| KRAS G12D/G12V/G12R targets + spectrum | Task 3 |
| Both-strand protospacer enumeration, PAM panel, SNV position, mutant-specific PAM | Task 4 |
| Discrimination margin (primary), CFD-style model, mechanism | Task 5 |
| On-target proxy (tie-break, binding caveat) | Task 5 |
| Genome-wide off-target search, specificity, veto, WT-allele as primary off-target | Tasks 5–6, 8 |
| Transparent lexicographic ranking + undesignable flag | Task 7 |
| Patient-coverage math | Task 7 |
| Reproducible run, cached/optional heavy step, README | Tasks 8–9 |
| ≤2-page writeup, trade-offs, honest limitations | Task 9 |
| Tests for enumeration / scoring / coverage | Tasks 4–7 |

Gaps: none. Note the spec's "sensitivity analysis" is realized as the `--margin-threshold` knob plus the transparent qualify-then-sort ranking (re-running at different thresholds shows stability); a dedicated sensitivity sweep is called out in the writeup as a quick extension rather than a separate module, to keep scope to one plan.

**2. Placeholder scan:** The only bracketed `<...>` values live in `writeup/writeup.md` and are *intentional* — Task 9 Step 2 instructs filling them from real results. No `TODO`/`TBD`/"implement later" in code tasks; every code step shows complete code.

**3. Type consistency:** `Guide`, `Target`, `DiscriminationResult`, `RankedGuide`, `OffTarget` field names are used identically across tasks. `evaluate_guide`/`rank_guides`/`recommend_per_allele`/`coverage_curve`/`run` signatures match their call sites in Tasks 7–8. Position convention (index 0 = PAM-distal) is consistent across `enumerate.py`, `cfd.py`, and `score.py`.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-06-24-pdac-crispr-guide-design.md`.
