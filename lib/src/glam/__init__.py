# Imports ======================================================================

# Standard Library
from pathlib import Path
from typing import Pattern

# Dependencies
import pyteomics.fasta

# Private Modules
import glam._lib as _lib

# Constants ====================================================================

GLYCOSYLATION_MOTIFS: dict[str, str] = {"N": r"N[^P][TS]"}

# Functions ====================================================================


# FIXME: Worth looking closely at argument names here, and make sure to update
# the docstring!
def generate_glycopeptides(
    fasta: str,
    enzyme: str | Pattern[str],
    motif: str | Pattern[str],
    glycan_database: str,
    missed_cleavages: int = 0,
    min_length: int | None = None,
    max_length: int | None = None,
    semi_enzymatic: bool = False,
) -> list[tuple[str, str]]:
    """Generates peptides from a protein sequence using the specified rule.

    Parameters
    ----------
    seq: str,
        The protein sequence to digest — a `str` of single-letter amino acids
    rule: str | Pattern[str], default="trypsin",
        The "rule" / enyzme used to calculate cleavage sites. A number of
        enzymes are built-in (see pyteomics.parser.expasy_rules and
        pyteomics.parser.psims_rules), but a custom (optionally pre-compiled)
        regex can also be supplied
    missed_cleavages: int, default=0,
        The maximum number of missed cleavages to allow during digestion
    min_length: int | None, default=None,
        The minimum length peptide to include in the digest output
    max_length: int | None, default=None,
        The maximum length peptide to include in the digest output
    semi: bool, default=False,
    """

    # entries = pyteomics.fasta.read(fasta)
    # glycans = _lib.load_glycans(glycan_database)

    # def generate(entry: FastaEntry) -> tuple[str,str]:
    #     seq = entry.sequence
    #     pass

    # # test_glycans = Path("./data/chlamy_glycans.csv").read_text()

    # # # Digest the protein + get candidate peptides
    # # peptides = _lib.filter_glycopeptides(_lib.digest_protein(spike_protein))

    # # # Load glycans
    # # glycans = _lib.load_glycans(test_glycans)

    # # # Generate glycopeptides
    # # glycopeptides = _lib.generate_glycopeptides(peptides, glycans)

    # return [generate(entry) for entry in entries]
    return []
