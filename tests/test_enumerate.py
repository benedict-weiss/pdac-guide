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
