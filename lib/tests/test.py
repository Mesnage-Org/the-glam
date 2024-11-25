from pathlib import Path

# Test Data  ===================================================================

SPIKE_PROTEIN: str = Path("tests/data/algal_spike.faa").read_text()
CHLAMY_GLYCANS: str = Path("tests/data/chlamy_glycans.csv").read_text()

# Integration Tests ============================================================


# FIXME: Write regression tests using pytest-regtest
def test_generate_glycopeptides_fasta() -> None:
    pass
