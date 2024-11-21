# Imports =====================================================================

import re
from typing import Pattern
from pyteomics.parser import cleave

# Constants ===================================================================

GLYC_MOTIFS: dict[str, str] = {"N": "N[^P][TS]"}


# Functions ===================================================================


def digest_protein(
    seq: str,
    rule: str | Pattern[str] = "trypsin",
    missed_cleavages: int = 0,
    min_length: int | None = None,
    max_length: int | None = None,
    semi: bool = False,
) -> set[str]:
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

    return cleave(seq, rule, missed_cleavages, min_length, max_length, semi)


def filter_glycopeptides(
    peptides: set[str], glycosylation_motif: str | Pattern[str] = GLYC_MOTIFS["N"]
) -> set[str]:
    return {p for p in peptides if re.search(glycosylation_motif, p)}
