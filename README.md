# PDAC CRISPR Guide Design (KRAS codon 12)

Mutation-specific guide-RNA design for a conditional CRISPR therapy: recognition
of a cancer mutation triggers a kill switch, so the guide carries both efficacy
and safety. This repo designs and ranks guides for the KRAS codon-12 alleles
(G12D / G12V / G12R) that dominate PDAC.

## Quickstart

```bash
pip install -e ".[dev]"
pytest                        # run the test suite
python3 scripts/run_analysis.py
```

Use `python3` (the script needs Python 3.11+). Outputs: `results/recommendations.tsv`
and `results/figures/*.png`. Without a reference the off-target scan is skipped and
the `specificity` column stays empty.

## Genome-wide off-target search (optional, heavier)

The core run scores discrimination and the wild-type-allele safety check without
any download. For a genome-scale off-target scan:

```bash
python3 scripts/fetch_data.py --out data/chr12.fa     # downloads chr12 (~129 MB, contains KRAS)
python3 scripts/run_analysis.py --reference data/chr12.fa
```

This fills in the `specificity` column and applies the safety veto. The scan is
vectorized (NumPy), so chr12 takes roughly 100 seconds; a progress bar shows on a
terminal. Chromosome 12 is a documented scope limit; full-genome (all chromosomes,
or Cas-OFFinder if installed) is the extension.

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

The CFD model is a documented approximation; on-target predictors are
cutting-trained while the trigger is binding; the off-target scan covers chr12
only (full-genome plus population variants are the extension); and locus
accessibility, WT-KRAS LOH, delivery, and wet-lab validation are out of scope.
