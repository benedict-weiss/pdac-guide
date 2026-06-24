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
G12R (~17%) of KRAS-mutant tumors (COSMIC/TCGA-PAAD/GENIE); each is a distinct
DNA sequence needing its own guide.

## Method: discrimination first
For each allele we enumerate every candidate protospacer (both strands, NGG PAM,
extensible to other PAMs) whose recognition overlaps the mutation, then score:

1. **Discrimination margin** Δ = score(mutant) − score(WT), using a CFD-style
   position-weighted mismatch model (Doench et al. 2016). A mismatch in the
   PAM-proximal seed, or a mutation that gives the mutant a PAM the WT lacks,
   drives WT recognition to ~0.
2. **Off-target safety veto** — a genome/chr12 scan, CFD-aggregated, vetoing any
   guide with a high-CFD off-target; the WT allele is checked explicitly. *(Not
   run in the default configuration — requires `--reference`; see Limitations.)*
3. **On-target proxy** — a documented guardrail (GC, homopolymers, Pol III
   terminator), used only to break ties. It is *not* a validated predictor.

Ranking is lexicographic and inspectable, not a black box.

## Recommendation
| Allele | Spacer (5'→3') | Strand | Mechanism | Margin | Coverage added |
|--------|----------------|--------|-----------|--------|----------------|
| G12D | `CTTGTGGTAGTTGGAGCTGA` | `+` | `seed_mismatch` | 0.98 | ~36% |
| G12V | `CTTGTGGTAGTTGGAGCTGT` | `+` | `seed_mismatch` | 0.98 | ~27% |
| G12R | `CTTGTGGTAGTTGGAGCTCG` | `+` | `seed_mismatch` | 0.97 | ~15.3% |

All three guides share a downstream `TGG` PAM (GRCh38 KRAS locus, chr12). The
codon-12 SNV lands in the PAM-proximal seed (positions 1–8 from PAM), where
mismatches most strongly reduce Cas9 binding. This places WT recognition at
0.02–0.03 (near-zero) for all three guides.

The three-guide set addresses ~78.3% of PDAC patients (frequency × 0.90
probability-of-recognition, cumulative). Marginal coverage per allele: 36% (G12D) +
27% (G12V) + 15.3% (G12R). No allele in this set was flagged undesignable.

**Important:** discrimination (mutant vs. WT) is strong for all three alleles.
Genome-wide off-target safety has **not** been evaluated in the default run — the
`specificity` column is empty because no reference genome was provided. Evaluating
it is a two-step flow — `scripts/fetch_data.py --out data/chr12.fa` to download
chr12, then `scripts/run_analysis.py --reference data/chr12.fa` to re-run the
pipeline with the scan. That scan is the immediate next step and an explicit
limitation of this report.

## Trade-offs I weighed
- **Discrimination vs coverage:** chasing rarer alleles adds patients but often
  worse margins; I prioritized clean discrimination over coverage.
- **Seed mismatch vs mutant-specific PAM:** PAM-differential guides give the
  cleanest separation but are not always available per allele. All three top
  guides here exploit seed mismatch with a shared `TGG` PAM.
- **Reproducibility vs scale:** the default run needs no download; genome-scale
  off-target is an explicit, cached extension.

## Honest limitations / with more time
- **Model:** CFD (Doench et al. 2016) is an approximation of *cutting* data;
  the trigger is *binding*. I would substitute the full empirical CFD/CRISPRoff-
  style tables and a binding-specific model, and validate in isogenic
  mutant-vs-WT reporter lines.
- **Genome-wide off-targets not evaluated:** the default run skips the chr12
  scan (no `--reference`). True safety requires all chromosomes plus common
  population variants — a guide can be specific against the reference but not
  against a patient whose WT allele carries a SNP.
- **WT-KRAS LOH / zygosity:** if a tumor loses the WT allele, discrimination
  matters less; copy-number context should feed target choice.
- **Locus accessibility & expression**, guide secondary structure, delivery/PK,
  and immunogenicity are out of scope here.
- **Off-target scope:** chr12 by default; full-genome plus population-variant
  sweeps are the extension, as documented in `scripts/fetch_data.py`.
