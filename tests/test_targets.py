from guide_design.targets import KRAS_REF_CONTEXT, Target, load_targets, TARGETS


def test_reference_context_codon12_is_ggt():
    # Codon 12 sits at 0-based indices 33-35 of the CDS context.
    assert len(KRAS_REF_CONTEXT) == 90
    assert KRAS_REF_CONTEXT[33:36] == "GGT"


def test_targets_loaded_for_three_alleles():
    alleles = {t.allele for t in TARGETS}
    assert alleles == {"G12D", "G12V", "G12R"}


def test_ref_base_matches_context():
    for t in TARGETS:
        assert KRAS_REF_CONTEXT[t.snv_index] == t.ref_base


def test_mut_seq_applies_single_substitution():
    g12d = next(t for t in TARGETS if t.allele == "G12D")
    mut = g12d.mut_seq
    assert mut[34] == "A"                 # G>A applied
    assert mut[:34] == KRAS_REF_CONTEXT[:34]   # nothing else changed
    assert mut[35:] == KRAS_REF_CONTEXT[35:]
    assert mut[33:36] == "GAT"            # codon 12 now Asp


def test_g12r_changes_first_codon_base():
    g12r = next(t for t in TARGETS if t.allele == "G12R")
    assert g12r.snv_index == 33
    assert g12r.mut_seq[33:36] == "CGT"   # Arg
