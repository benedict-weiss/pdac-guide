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
