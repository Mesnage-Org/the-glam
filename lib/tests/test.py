from pathlib import Path

from glam import GLYCOSYLATION_MOTIFS, generate_glycopeptides

# Test Data  ===================================================================

SPIKE_PROTEIN: str = Path("tests/data/algal_spike.faa").read_text()
CHLAMY_GLYCANS: str = Path("tests/data/chlamy_glycans.csv").read_text()

# Integration Tests ============================================================


def test_generate_glycopeptides_fasta() -> None:
    pass
