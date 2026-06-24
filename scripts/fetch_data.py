"""OPTIONAL: build a reference for genome-wide off-target search.

Downloads a reference FASTA (default: human chromosome 12 from Ensembl, which
contains KRAS) so off-target search runs against real sequence. Full-genome is
the documented extension (concatenate all chromosomes, or use Cas-OFFinder).

Usage:
    python scripts/fetch_data.py --out data/chr12.fa
    python scripts/run_analysis.py --reference data/chr12.fa

NOTE: chromosome 12 alone is a deliberate, documented scope limit -- true
off-targets can occur on any chromosome. State this in the writeup.

We pull from the Ensembl FTP bulk download (a gzipped whole-chromosome FASTA),
NOT the REST /sequence endpoint -- REST caps requests at 10 Mb, and chr12 is
~133 Mb, so REST returns HTTP 400 for a whole chromosome.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import urllib.request

# Ensembl FTP: gzipped soft-masked chromosome 12 FASTA (GRCh38). ~40 MB download.
ENSEMBL_CHR12 = (
    "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz"
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/chr12.fa")
    ap.add_argument("--url", default=ENSEMBL_CHR12)
    args = ap.parse_args()
    print(f"Downloading {args.url} -> {args.out} (this is the heavy, optional step)...")
    gz_path = args.out + ".gz"
    urllib.request.urlretrieve(args.url, gz_path)
    print(f"Decompressing {gz_path} -> {args.out} ...")
    with gzip.open(gz_path, "rb") as src, open(args.out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    print("Done. Re-run: python scripts/run_analysis.py --reference", args.out)


if __name__ == "__main__":
    main()
