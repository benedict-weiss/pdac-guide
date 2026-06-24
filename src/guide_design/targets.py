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
