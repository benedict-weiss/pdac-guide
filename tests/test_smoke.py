import guide_design


def test_package_imports_with_version():
    assert guide_design.__version__ == "0.1.0"
