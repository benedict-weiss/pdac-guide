# Designing a mutation-specific CRISPR therapy for PDAC

## What we're actually designing
FinalDose's molecule fires its payload when the guide RNA recognizes its target, which means the guide RNA basically *is* the therapy. Both the efficacy and the safety of the whole thing come down to that one sequence. So "guide design" is really three questions tangled together: what to go after, how to tell a mutant allele apart from a wild-type (WT) allele that differs by a single base, and how to keep it from firing somewhere else in the genome. And the off-target we should worry about most isn't some random genomic site. It's the WT allele sitting in every healthy cell in the patient.

## Why I targeted KRAS codon 12
About 90% of PDAC tumors carry a KRAS mutation, and those mutations pile up at one hotspot: codon 12. KRAS is an essential driver, so the tumor can't just delete it to dodge the therapy. That's a real advantage over TP53, CDKN2A, or SMAD4, which tend to get lost by deletion or show up only in subclones, both of which make them lousy point-mutation triggers. So I designed at codon 12 and made the case for the rest in prose.

Within codon 12, three alleles dominate: G12D (~40%), G12V (~30%), and G12R (~17%) of KRAS-mutant tumors (per COSMIC, TCGA-PAAD, and GENIE). Each one is a different DNA sequence, so each needs its own guide.

## How the method works: discrimination first
For each allele I enumerate every candidate protospacer whose recognition window overlaps the mutation (both strands, NGG PAM, and the code extends to other PAMs), then score it on three things:

1. **Discrimination margin**, Δ = score(mutant) − score(WT). I'm using a CFD-style position-weighted mismatch model (Doench et al. 2016). A mismatch in the PAM-proximal seed, or a mutation that hands the mutant a PAM the WT doesn't have, pushes WT recognition down toward zero.
2. **Off-target safety veto**: a genome (or chr12) scan, CFD-aggregated, that vetoes any guide with a high-CFD off-target. The WT allele gets checked explicitly. Worth flagging up front: this doesn't run in the default config. It needs `--reference`. See Limitations.
3. **On-target proxy**: a documented guardrail (GC content, homopolymers, Pol III terminator) that I only use to break ties. It is *not* a validated predictor, and I don't treat it as one.

The ranking is lexicographic and you can read it top to bottom. No black box.

## What I'd recommend
| Allele | Spacer (5'→3') | Strand | Mechanism | Margin | Specificity (chr12) | Coverage added |
|--------|----------------|--------|-----------|--------|---------------------|----------------|
| G12D | `CTTGTGGTAGTTGGAGCTGA` | `+` | `seed_mismatch` | 0.98 | 76.2 | ~36% |
| G12V | `CTTGTGGTAGTTGGAGCTGT` | `+` | `seed_mismatch` | 0.98 | 76.8 | ~27% |
| G12R | `CTTGTGGTAGTTGGAGCTCG` | `+` | `seed_mismatch` | 0.97 | 96.8 | ~15.3% |

All three guides share a downstream `TGG` PAM at the GRCh38 KRAS locus on chr12. The codon-12 SNV falls in the PAM-proximal seed (positions 1 to 8 from the PAM), which is exactly where mismatches hurt Cas9 binding the most. That puts WT recognition at 0.02 to 0.03, near zero, for all three.

Together the three guides cover roughly 78.3% of PDAC patients (allele frequency × the 0.90 KRAS-mutant fraction, summed up). Per allele that breaks down as 36% (G12D) + 27% (G12V) + 15.3% (G12R). None of the three came back undesignable.

For this report I ran the chr12 off-target scan (`scripts/fetch_data.py --out data/chr12.fa`, then `scripts/run_analysis.py --reference data/chr12.fa`). All three guides clear the safety veto: no single off-target on chr12 has a CFD above the 0.5 threshold, so none of them gets vetoed. The specificity column above is the MIT-style aggregate (100 means perfectly specific). G12R is the cleanest at 96.8. G12D and G12V sit lower, around 76, which means they pick up a heavier but still sub-veto off-target load on chr12.

The thing I won't oversell: this is chr12 only. A real safety sign-off needs every chromosome plus common population variants, so treat these numbers as a strong first pass, not a clearance. More on that in the limitations.

## Trade-offs I had to make
- **Discrimination vs. coverage.** Going after rarer alleles picks up more patients but usually at worse margins, so I chose clean discrimination over chasing coverage.
- **Seed mismatch vs. mutant-specific PAM.** PAM-differential guides give you the cleanest separation, but you can't always find one per allele. All three top guides here lean on a seed mismatch with a shared `TGG` PAM instead.
- **Reproducibility vs. scale.** The default run needs nothing downloaded. Genome-scale off-target search is a deliberate, cached extension rather than something baked into every run.

## Where this falls short, and what I'd do with more time
- **The model.** CFD (Doench et al. 2016) is fit to *cutting* data, but the trigger here is *binding*. I'd swap in the full empirical CFD / CRISPRoff-style tables and a binding-specific model, then validate in isogenic mutant-vs-WT reporter lines.
- **Off-target scan covers chr12 only.** I ran the chr12 scan for this report and all three guides passed, but real safety means all chromosomes plus common population variants. A guide can look specific against the reference and still misfire in a patient whose WT allele happens to carry a SNP. The default run (no `--reference`) skips the scan entirely.
- **WT-KRAS loss of heterozygosity / zygosity.** If a tumor has already lost its WT allele, discrimination matters less, so copy-number context really ought to feed back into target choice.
- **Out of scope here:** locus accessibility and expression, guide secondary structure, delivery and PK, and immunogenicity.
- **Off-target scope.** Chr12 by default; full-genome plus population-variant sweeps are the extension, documented in `scripts/fetch_data.py`.
