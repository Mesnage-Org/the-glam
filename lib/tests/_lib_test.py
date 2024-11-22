from textwrap import dedent
from pathlib import Path

import pytest
from glam import GLYCOSYLATION_MOTIFS
from glam._lib import (
    glycan_mass,
    load_glycans,
    digest_protein,
    filter_glycopeptides,
    peptide_mass,
    convert_to_csv,
)
from pyteomics import fasta

# Test Data  ===================================================================

SPIKE_PROTEIN: str = next(fasta.read("tests/data/algal_spike.faa")).sequence
GLYCANS: str = Path("tests/data/chlamy_glycans.csv").read_text()
GLYCANS_AND_MASSES: str = Path("tests/data/chlamy_glycans_and_masses.csv").read_text()

# Helper Functions ============================================================


def detrim(s: str) -> str:
    return dedent(s.strip("\r\n"))


# Unit Tests ===================================================================

def test_glycan_mass() -> None:
    iupac_structure = "Neu5Ac(a2-3)Gal6S(b1-3)[Neu5Ac(a2-6)]GalNAc"
    assert glycan_mass(iupac_structure) == 1045.2903546
    # byonic_composition = "Hex(4)HexNAc(2)Me(3)Hex(3)Pent(2)dHex(1)"
    # assert glycan_mass(byonic_composition) == 2010.72845026293

def test_load_glycans_with_masses() -> None:
    glycans = load_glycans(GLYCANS_AND_MASSES)
    assert len(glycans) == 52
    assert all(type(g) is str and type(m) is float for g, m in glycans)


def test_load_glycans_without_masses() -> None:
    glycans = load_glycans(GLYCANS)
    assert len(glycans) == 52
    assert all(type(g) is str and type(m) is float for g, m in glycans)
    expected_glycans = load_glycans(GLYCANS_AND_MASSES)
    things = [(abs(e[1] - o[1]), o[1], o[0]) for e, o in zip(sorted(list(expected_glycans)), sorted(list(glycans)))]
    for thing in sorted(things):
        print(thing)
    breakpoint()


def test_load_glycans_raises() -> None:
    no_header = detrim("""
    Hex(6)HexNAc(2)MeHex(0)Pent(0)dHex(1),1524.5335944644
    Hex(4)HexNAc(2)MeHex(2)Pent(2)dHex(0),1670.5915032738
    Hex(5)HexNAc(2)MeHex(1)Pent(2)dHex(0),1656.5758532096
    """)
    with pytest.raises(ValueError) as e:
        load_glycans(no_header)
    assert (
        str(e.value)
        == "The glycan file must contain the columns: ['Glycan', 'Monoisotopic Mass']"
    )


def test_digest_protein() -> None:
    expected_tryptic_peptides = {
        "AAEIR",
        "AGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTK",
        "AHFPR",
        "ALTGIAVEQDK",
        "ASANLAATK",
        "CTLK",
        "CVNFNFNGLTGTGVLTESNK",
        "CYGVSPTK",
        "DFGGFNFSQILPDPSKPSK",
        "DIADTTDAVR",
        "DISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYR",
        "DLICAQK",
        "DLPQGFSALEPLVDLPIGINITR",
        "DPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWR",
        "EELDK",
        "EFVFK",
        "EGVFVSNGTHWFVTQR",
        "EIDR",
        "FASVYAWNR",
        "FDNPVLPFNDGVYFASTEK",
        "FLPFQQFGR",
        "FNGIGVTQNVLYENQK",
        "FNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYR",
        "FPNITNLCPFGEVFNATR",
        "FQTLLALHR",
        "GDEVR",
        "GIYQTSNFR",
        "GVYYPDK",
        "GWIFGTTLDSK",
        "GYHLMSFPQSAPHGVVFLHVTYVPAQEK",
        "HTPINLVR",
        "IADYNYK",
        "IQDSLSSTASALGK",
        "ISNCVADYSVLYNSASFSTFK",
        "IYSK",
        "K",
        "LDPPEAEVQIDR",
        "LFR",
        "LIANQFNSAIGK",
        "LITGR",
        "LNDLCFTNVYADSFVIR",
        "LNEVAK",
        "LPDDFTGCVIAWNSNNLDSK",
        "LQDVVNQNAQALNTLVK",
        "LQSLQTYVTQQLIR",
        "MAR",
        "MNLTTR",
        "MSECVLGQSK",
        "NFTTAPAICHDGK",
        "NFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFK",
        "NHTSPDVDLGDISGINASVVNIQK",
        "NIDGYFK",
        "NK",
        "NLNESLIDLQELGK",
        "NLR",
        "NNK",
        "NTQEVFAQVK",
        "QGNFK",
        "QIAPGQTGK",
        "QIYK",
        "QLSSNFGAISSVLNDILSR",
        "QYGDCLGDIAAR",
        "R",
        "SFIEDLLFNK",
        "SFTVEK",
        "SNIIR",
        "SNLKPFER",
        "SSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTK",
        "STNLVK",
        "SWMESEFR",
        "SYLTPGDSSSGWTAGAAAYYVGYLQPR",
        "TFLLK",
        "TGALLLVALALAGCAQACR",
        "TPPIK",
        "TQLPPAYTNSFTR",
        "TQSLLIVNNATNVVIK",
        "TSVDCTMYICGDSTECSNLLLQYGSFCTQLNR",
        "VCEFQFCNDPFLGVYYHK",
        "VDFCGK",
        "VFR",
        "VGGNYNYLYR",
        "VQPTESIVR",
        "VTLADAGFIK",
        "VVVLSFELLHAPATVCGPK",
        "VYSSANNCTFEYVSQPFLMDLEGK",
        "VYSTGSNVFQTR",
        "WGSHHHHHHHHSPSPSPSPSPSPSPSPSPSPSPSPSPSPSPSPSPSPSPSPSGYPYDVPDYA",
        "YEQYIK",
        "YFK",
        "YNENGTITDAVDCALDPLSETK",
    }
    tryptic_peptides = digest_protein(SPIKE_PROTEIN)
    assert len(tryptic_peptides) == 90
    assert tryptic_peptides == expected_tryptic_peptides


def test_filter_glycopeptides() -> None:
    expected_glycopeptides = {
        "AGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTK",
        "DFGGFNFSQILPDPSKPSK",
        "DLPQGFSALEPLVDLPIGINITR",
        "DPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWR",
        "EGVFVSNGTHWFVTQR",
        "FPNITNLCPFGEVFNATR",
        "MNLTTR",
        "NFTTAPAICHDGK",
        "NFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFK",
        "NHTSPDVDLGDISGINASVVNIQK",
        "NLNESLIDLQELGK",
        "SSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTK",
        "TQSLLIVNNATNVVIK",
        "VYSSANNCTFEYVSQPFLMDLEGK",
        "YNENGTITDAVDCALDPLSETK",
    }
    tryptic_peptides = digest_protein(SPIKE_PROTEIN)
    glycopeptides = filter_glycopeptides(tryptic_peptides, GLYCOSYLATION_MOTIFS["N"])
    assert len(glycopeptides) == 15
    assert glycopeptides == expected_glycopeptides


def test_peptide_mass() -> None:
    mass = peptide_mass("PEPTIDE")
    assert mass == 799.3599640267099


def test_peptide_mass_raises() -> None:
    with pytest.raises(ValueError) as e:
        peptide_mass("PEPTIXE")
    assert (
        str(e.value)
        == "Unknown amino acid residue found in 'PEPTIXE': No mass data for residue: X"
    )


def test_build_glycopeptides() -> None:
    # FIXME: Check this is working with Tia!
    pass


def test_convert_to_csv() -> None:
    glycopeptides = {
        ("A", 42.123456789),
        ("B", 128.123456789),
        ("C", 1337.123456789),
    }
    expected_csv = detrim("""
        Structure,Monoisotopicmass
        A,42.123457
        B,128.123457
        C,1337.123457
    """).replace("\n", "\r\n")
    csv = convert_to_csv(glycopeptides)
    print(csv)
    assert csv == expected_csv
