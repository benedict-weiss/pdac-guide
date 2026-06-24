# PDAC CRISPR Guide-RNA Design — Design Spec

- **Date:** 2026-06-24
- **Status:** Approved (design); pending implementation plan
- **Brief:** `brief.md` — WP3 take-home, "CRISPR guide design for a cancer indication" (PDAC)
- **Time-box:** ~1–2 days

---

## 1. Problem framing

FinalDose's therapy is a **mutation-specific conditional CRISPR** system: a guide RNA directs a
molecule to a DNA sequence, and **recognition itself is the trigger** that releases a cell-killing
payload. Cancer cells carrying the targeted mutation die; healthy cells carrying the wild-type (WT)
sequence are left untouched.

The consequence the whole project hinges on: **both efficacy and safety live in the guide.**
Wherever recognition happens is where the payload fires — on the intended mutation (therapy) or
anywhere it shouldn't (a dead healthy cell). So "guide design" is really three coupled problems:

1. **Target selection** — which somatic mutation marks PDAC cells well enough to aim at?
2. **Discrimination** — most PDAC drivers are single-nucleotide variants (SNVs). The hardest and most
   safety-critical "off-target" is the **WT allele present in every healthy cell, one base away.**
3. **Genome-wide safety** — because binding fires the payload, any sufficiently similar genomic site
   is a potential healthy-cell kill, so off-target enumeration + scoring is mandatory, not optional.

The deliverable is a **defensible decision path**: "treat PDAC" → pick target(s) → design candidate
guides → score them (discrimination, off-target safety, on-target efficacy) → judge whether the
design is actually good → recommend, with honest limitations.

## 2. Goals / non-goals

**Goals**
- A reproducible, real-but-lightweight pipeline that designs and ranks guides for the **KRAS
  codon-12 allelic series (G12D, G12V, G12R)**.
- A transparent, inspectable scoring model with discrimination as the primary objective.
- A patient-coverage analysis tying the guide set to the fraction of PDAC patients addressed.
- A ≤2-page writeup carrying the reasoning and honest limitations.

**Non-goals**
- Designing guides for non-KRAS drivers (TP53/CDKN2A/SMAD4) — argued in prose during target
  selection, **not** coded. They are frequently subclonal or copy-number events, a poor fit for a
  point-mutation trigger.
- A tuned black-box predictor. The brief explicitly prefers honest limitations over a tuned model.
- Wet-lab validation (scoped as future work).
- A production web service / UI.

## 3. Locked assumptions

1. **Mechanism = binding / R-loop-triggered, enzyme-agnostic.** Concrete anchor is SpCas9 geometry
   (20-nt spacer, NGG PAM), but **PAM is treated as a design *variable*** — we also evaluate
   alternative PAMs, because a mutation that *creates or destroys* a PAM is the single strongest
   discrimination lever. On-target *cutting*-efficiency predictors are used only as a soft efficacy
   proxy, with that caveat stated explicitly.
2. **Reference = GRCh38, canonical KRAS transcript**; mutation spectrum from public PDAC data
   (COSMIC / TCGA / GENIE), curated into a small cited table.
3. **The WT allele is the primary safety off-target**, evaluated explicitly alongside the
   genome-wide search.

## 4. Scientific approach (decision tree)

1. **Why KRAS.** Truncal/clonal in ~90% of PDAC; a single recurrent hotspot codon (12) makes it
   tractable; it is an essential oncogenic driver, so escape by target loss is unlikely (unlike a
   passenger). Contrast with TP53 (often LOH/deletion), CDKN2A (deletion), SMAD4 (deletion/subclonal).
2. **Which alleles.** KRAS codon-12 spectrum in PDAC, approx. as fraction of KRAS-mutant tumors:
   G12D ≈ 40%, G12V ≈ 30%, G12R ≈ 17%, remainder (G12C, Q61, etc.) smaller. Exact figures pinned in
   the data table. We design one guide per allele because **each allele is a different DNA sequence**.
   - Codon 12 reference = GGT (Gly). G12D = GGT→GAT (mid base), G12V = GGT→GTT (mid base),
     G12R = GGT→CGT (first base). The changed base sits at different positions, which changes where
     it can land relative to candidate PAMs/seed — the crux of the per-allele comparison.
3. **Design engine (Approach A, discrimination-first):** enumerate protospacers × PAM panel × both
   strands in a window around each variant → score discrimination margin Δ(mutant − WT) **first** →
   apply genome-wide off-target safety **veto** → rank survivors by on-target proxy → recommend one
   guide per allele.
4. **Coverage:** combine per-allele recommendations into a coverage curve over patients.
5. **Judge it:** per-allele achievable margin, off-target worst case, WT-allele safety, sensitivity
   analysis, and honest limitations.

## 5. Architecture

```
pdac-guide/
├── README.md                  # reproducible run + decision tree + limitations
├── pyproject.toml
├── data/
│   ├── kras_locus.fa          # small GRCh38 slice around KRAS codon 12 (bundled)
│   ├── pdac_kras_spectrum.tsv # curated allele freqs, cited (COSMIC/TCGA/GENIE)
│   └── offtarget_cache/       # cached genome-wide hits → reviewer needs no full-genome download
├── src/guide_design/
│   ├── targets.py     # target defs (gene, GRCh38 coords, ref/alt, freq); build mut vs WT context
│   ├── enumerate.py   # protospacers: both strands, window, PAM panel; SNV position; PAM create/destroy flags
│   ├── cfd.py         # vendored CFD mismatch/PAM weight tables (cited)
│   ├── score.py       # discrimination score + on-target proxy
│   ├── offtarget.py   # genome-wide search (Cas-OFFinder if present, else bundled k-mismatch) + CFD aggregation
│   ├── rank.py        # combine → ranked table + recommended guide/allele + rationale
│   └── coverage.py    # patient-coverage math from allele freqs
├── scripts/
│   ├── fetch_data.py      # regenerate bundled data from public sources (optional, documented)
│   └── run_analysis.py    # end-to-end → results/ + figures
├── results/               # recommendations.tsv + figures/*.png
├── tests/                 # enumerate / score / coverage invariants
└── writeup/writeup.md     # ≤2 pages
```

**Reproducibility contract:** `pip install -e . && python scripts/run_analysis.py` runs end-to-end
in minutes using cached off-target data — **no multi-GB genome download required**. `fetch_data.py`
documents and regenerates the bundled data from public sources; the full-genome off-target path is
documented for anyone who wants to rerun it from scratch.

## 6. Module specifications

### `targets.py`
- Defines each target: gene, GRCh38 coordinates, strand, reference/alt base(s), HGVS, allele
  frequency (from the spectrum table).
- Builds the local **mutant** and **WT** sequence contexts around each variant from `kras_locus.fa`.
- Exact coordinates are pinned here / in `data` (not hard-coded into prose) to avoid drift.

### `enumerate.py`
- Given a variant + flanking context, enumerate candidate protospacers on **both strands** across a
  window (default ±~25 nt) for a **PAM panel** (default: SpCas9 `NGG`; configurable to add e.g. `NG`,
  Cas12a `TTTV`).
- For each candidate record: spacer sequence, strand, PAM, cut/recognition frame, and the **position
  of the SNV within the protospacer** (PAM-proximal = seed).
- Detect when the SNV makes the **PAM differ between mutant and WT**. The strongest discrimination
  class is a **mutant-specific PAM** (PAM present in the mutant context, absent in WT → WT cannot be
  bound at all). The reverse (PAM present in WT but destroyed in the mutant) is *unusable* for
  targeting the mutant and is flagged/excluded, not rewarded.
- Only candidates whose protospacer actually overlaps the variant base are retained (a guide that
  doesn't cover the SNV can't discriminate).

### `cfd.py`
- Vendored, cited CFD (Cutting Frequency Determination) position- and identity-weighted mismatch
  tables, plus PAM-activity weights. Pure data + a lookup function; no network dependency.

### `score.py`
- **`discrimination(guide)`** — primary. Compute a CFD-style recognition score of the guide against
  (a) the mutant context and (b) the WT context. Return:
  - `score_mut`, `score_wt`, `margin = score_mut − score_wt`, and the **mechanism**
    (`seed_mismatch` vs `pam_differential` — a mutant-specific PAM drives `score_wt ≈ 0`).
  - A guide passes the discrimination gate iff `margin ≥ τ` (default threshold configurable).
- **`on_target_proxy(guide)`** — secondary, tie-break only. A documented activity heuristic
  (GC content, homopolymer runs, position features; Rule Set 2 / Azimuth if a clean library is
  available). **Caveat carried in output:** trained on cutting, not binding.

### `offtarget.py`
- Genome-wide enumeration of sites within ≤ **N=4** mismatches (default).
  - Use **Cas-OFFinder** if installed; otherwise a **bundled pure-Python k-mismatch search** over
    GRCh38, with results **cached** in `data/offtarget_cache/`.
- Aggregate hits into an MIT-style specificity score via CFD.
- **Hard veto** any guide with a high-CFD off-target (configurable threshold), with extra weight on
  essential / healthy-critical loci.
- **Always evaluate the WT KRAS allele as an explicit site** — it is a 1-mismatch "off-target"
  present in every healthy cell and is the safety lynchpin.

### `rank.py`
- Lexicographic ranking: `margin ≥ τ` **AND** off-target veto pass → then sort by `on_target_proxy`.
- Emits a ranked table per allele, the recommended guide, and a human-readable rationale
  (why it won, which discrimination mechanism, off-target worst case).
- Includes a **sensitivity pass**: re-rank under perturbed thresholds/weights and report whether the
  top pick is stable. No single opaque composite score.

### `coverage.py`
- From the allele spectrum, compute the patient-coverage curve for guide sets
  {G12D} → {G12D, G12V} → {G12D, G12V, G12R}, plus the marginal gain per added guide.
- KRAS-mutant ≈ 90% of PDAC; the three-guide set reaches roughly ~75–78% of all PDAC patients
  (exact value from the data table).

## 7. Validation — "is the design actually good?"

- **Per-allele:** best achievable discrimination margin; explicitly flag any allele that is
  effectively **undesignable** (no PAM in range / SNV stuck PAM-distal) as an honest negative result
  rather than forcing a weak guide.
- **Specificity:** genome-wide off-target count + worst-case CFD + the explicit WT-allele check.
- **Coverage:** the curve + marginal value of each guide.
- **Sensitivity:** stability of the recommendation under threshold/weight perturbation.

## 8. Testing strategy

- `test_enumerate.py` — strand handling, PAM matching, SNV-position indexing, PAM create/destroy
  detection, rejection of protospacers that don't cover the SNV.
- `test_score.py` — discrimination monotonicity (seed mismatch ⇒ larger margin than PAM-distal),
  mutant-specific PAM ⇒ near-maximal margin (`score_wt ≈ 0`), CFD lookup edge cases.
- `test_coverage.py` — coverage math and marginal-gain arithmetic against hand-computed values.

## 9. Data sources & provenance

- **GRCh38** KRAS locus slice (bundled FASTA) — provenance + coordinates recorded in `fetch_data.py`.
- **PDAC KRAS allele spectrum** — curated from COSMIC / TCGA-PAAD / AACR GENIE, cited in the TSV and
  the writeup. Figures are approximate and the citation makes the source explicit.
- **CFD weights** — vendored from the published Doench et al. tables, cited.

## 10. Deliverables (mapping to the brief)

- **Code + short README** — reproducible one-command run, cached off-target data, documented
  full-genome path. (Brief: "Code — reproducible, short README.")
- **`writeup/writeup.md` (≤2 pages)** — approach, recommendation (the guide set), trade-offs, honest
  limitations, and the decision tree. (Brief: "Writeup ≤2 pages.")
- **Tests** — encode the engine's core invariants.

## 11. Limitations & future work (to be written honestly in the writeup)

- **Cutting-vs-binding model gap** — CFD and on-target predictors are trained on Cas9 *cutting*;
  our trigger is *binding/R-loop formation*. We use them as the best available proxy and say so.
- **Locus accessibility & expression** — chromatin state and KRAS expression at the target are not
  modeled; recognition in vivo depends on accessibility.
- **Tumor zygosity & LOH of WT KRAS** — if a tumor has lost the WT allele, discrimination matters
  less; if a healthy tissue context differs, it matters more. Not modeled here.
- **Guide RNA secondary structure**, delivery / PK, and immunogenicity — out of scope.
- **No wet-lab readout** — the real test is a reporter/kill assay in isogenic mutant-vs-WT lines;
  scoped as the next step.

## 12. Build sequence (feeds the implementation plan)

1. Scaffold (`pyproject.toml`, package, dirs) + bundled data (`targets`, KRAS locus, spectrum).
2. `enumerate.py` + tests.
3. `cfd.py` + `score.py` (discrimination) + tests.
4. `offtarget.py` + caching.
5. `rank.py` + `coverage.py` + tests.
6. `scripts/run_analysis.py` + figures.
7. `README.md` + `writeup/writeup.md`.
8. Verification pass: run end-to-end, sanity-check recommendations (and any literature-known KRAS
   allele-specific guides) for plausibility.

## 13. Defaults chosen (override points)

- Bundle a pure-Python off-target search; use Cas-OFFinder only if already installed.
- Off-target mismatch threshold **N = 4**.
- Vendor published CFD weight tables (cited) rather than depend on an external package.
- Enumeration window ±~25 nt around the variant; PAM panel defaults to SpCas9 `NGG`, configurable.
