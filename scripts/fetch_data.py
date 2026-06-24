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
