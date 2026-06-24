from guide_design.seq import revcomp, pam_matches


def test_revcomp_basic():
    assert revcomp("ATGC") == "GCAT"


def test_revcomp_is_involution():
    s = "ACGTACGTTTGGCCAA"
    assert revcomp(revcomp(s)) == s


def test_pam_matches_ngg():
    assert pam_matches("TGG") is True
    assert pam_matches("AGG") is True
    assert pam_matches("TAG") is False   # 2nd base not G


def test_pam_matches_length_guard():
    assert pam_matches("GG") is False     # wrong length


def test_pam_matches_case_insensitive():
    assert pam_matches("tgg") is True        # lowercase seq3 (default uppercase pattern)
    assert pam_matches("TGG", "ngg") is True  # lowercase pattern wildcard still matches
    assert pam_matches("TAG", "ngg") is False  # lowercase pattern, real mismatch
