# Imports =====================================================================

# Standard Library
from io import StringIO
from typing import Pattern
import csv
import itertools
import re

# Dependencies
import glycowork.motif.tokenization as glycowork
from pyteomics.auxiliary import PyteomicsError
import pyteomics.parser
import pyteomics.mass

# Constants ===================================================================

WATER_MASS: float = 18.01056468403

# Functions ===================================================================


def glycan_mass(glycan: str) -> float:
    try:
        return glycowork.glycan_to_mass(glycan)
    except KeyError:
        pass

    try:
        mass = glycowork.composition_to_mass(glycan)
        # NOTE: Many nonsense compositions just return the mass of water — if
        # the result was a water mass, then whatever we were given wasn't a
        # valid composition, so skip the return
        if mass != 18.0105546:
            return mass
    except IndexError:
        pass

    raise ValueError(f"Invalid glycan structure / composition: '{glycan}'")


def load_glycans(glycan_csv: str) -> set[tuple[str, float]]:
    expected_cols: list[str] = ["Glycan", "Monoisotopic Mass"]

    try:
        rows = csv.reader(StringIO(glycan_csv))
    except Exception as e:
        raise ValueError(f"Failed to read glycan CSV: {e}")

    header = next(rows)

    if header == expected_cols:
        return {(g, float(m)) for g, m in rows}
    elif header == expected_cols[:1]:
        return {(g, glycan_mass(g)) for (g,) in rows}
    else:
        raise ValueError(f"The glycan file must contain the columns: {expected_cols}")


def digest_protein(
    seq: str,
    rule: str | Pattern[str] = "trypsin",
    missed_cleavages: int = 0,
    min_length: int | None = None,
    max_length: int | None = None,
    semi: bool = False,
) -> set[str]:
    return pyteomics.parser.cleave(
        seq, rule, missed_cleavages, min_length, max_length, semi
    )


def filter_glycopeptides(
    peptides: set[str], glycosylation_motif: str | Pattern[str]
) -> set[str]:
    return {p for p in peptides if re.search(glycosylation_motif, p)}


def peptide_mass(peptide: str) -> float:
    try:
        return pyteomics.mass.fast_mass(peptide)
    except PyteomicsError as e:
        raise ValueError(
            f"Unknown amino acid residue found in '{peptide}': {e.message}"
        )


def build_glycopeptides(
    peptides: set[str], glycans: set[tuple[str, float]]
) -> set[tuple[str, float]]:
    def build(peptide: str, glycan: tuple[str, float]) -> tuple[str, float]:
        glycan_name, glycan_mass = glycan
        glycopeptide_sequence = glycan_name + peptide
        # This is a condensation reaction, so remember to take away ~18 Da
        glycopeptide_mass = glycan_mass + peptide_mass(peptide) - WATER_MASS
        return (
            glycopeptide_sequence,
            glycopeptide_mass,
        )

    return {build(p, g) for p, g in itertools.product(peptides, glycans)}


def convert_to_csv(glycopeptides: set[tuple[str, float]]) -> str:
    csv_str = StringIO()
    writer = csv.writer(csv_str)
    writer.writerow(["Structure", "Monoisotopicmass"])

    for name, mass in sorted(glycopeptides):
        mass = round(mass, 6)
        writer.writerow([name, mass])

    return csv_str.getvalue()
