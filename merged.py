# Imports =====================================================================

# Standard Library
from io import StringIO
from pathlib import Path
from typing import Pattern
import csv
import itertools
import re

# Dependencies
from pyteomics.auxiliary import PyteomicsError
from pyteomics import parser as pyt_parser, mass as pyt_mass, fasta as pyt_fasta
import pandas as pd

# Constants ===================================================================

GLYC_MOTIFS: dict[str, str] = {"N": "N[^P][TS]"}
WATER_MASS: float = 18.010565

# Functions ===================================================================


def load_glycans(glycan_csv: str) -> set[tuple[str, float]]:
    expected_cols: list[str] = ["Glycan", "Monoisotopic Mass"]

    try:
        df = pd.read_csv(StringIO(glycan_csv))
    except Exception as e:
        raise ValueError(f"Failed to read glycan CSV: {e}")

    if not all(col in expected_cols for col in df.columns):
        raise ValueError(f"The glycan file must contain the columns: {expected_cols}")

    return set(df.itertuples(index=False, name=None))


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

    return pyt_parser.cleave(seq, rule, missed_cleavages, min_length, max_length, semi)


def filter_glycopeptides(
    peptides: set[str], glycosylation_motif: str | Pattern[str] = GLYC_MOTIFS["N"]
) -> set[str]:
    return {p for p in peptides if re.search(glycosylation_motif, p)}


def peptide_mass(peptide: str) -> float:
    try:
        return pyt_mass.fast_mass(peptide)
    except PyteomicsError as e:
        raise (f"Unknown amino acid residue in '{peptide}': {e}")


# Function to generate glycopeptides and calculate their monoisotopic mass
# Taking into consideration when a glycan attaches to a peptide a covalent bond is formed.
# The bond releases H₂O so reducing the overall mass by 18.01528 Da

def generate_glycopeptides(
    peptides: set[str], glycans: set[tuple[str, float]]
) -> set[tuple[str, float]]:
    def build_glycopeptide(
        peptide: str, glycan: tuple[str, float]
    ) -> tuple[str, float]:
        glycan_name, glycan_mass = glycan
        glycopeptide_sequence = glycan_name + peptide
        glycopeptide_mass = glycan_mass + peptide_mass(peptide) - WATER_MASS
        return (
            glycopeptide_sequence,
            glycopeptide_mass,
        )

    return {build_glycopeptide(p, g) for p, g in itertools.product(peptides, glycans)}


# Main function
def main():
    spike_protein = pyt_fasta.read("data/algal_spike.faa").next().sequence
    test_glycans = Path("./data/chlamy_glycans.csv").read_text()

    # Digest the protein + get candidate peptides
    peptides = filter_glycopeptides(digest_protein(spike_protein))

    # Load glycans
    glycans = load_glycans(test_glycans)

    # Generate glycopeptides
    glycopeptides = generate_glycopeptides(peptides, glycans)

    # Write CSV to str
    with StringIO() as output:
        writer = csv.writer(output)
        for name, mass in glycopeptides:
            mass = round(mass, 6)
            writer.writerow([name, mass])
        print(output.getvalue())


if __name__ == "__main__":
    main()
