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
