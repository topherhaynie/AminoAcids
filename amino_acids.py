"""Class for managing amino acids."""

from enum import Enum
from typing import NamedTuple


class AminoAcidProperties(NamedTuple):
    """Properties of an amino acid."""

    full_name: str
    three_letter_code: str
    one_letter_code: str
    polar: bool
    group: str  # e.g., alkyl, aromatic, etc.
    charge: str  # e.g., neutral, positive, negative
    alt_name: str = ""

    @classmethod
    def to_dict(cls):
        """Return a dictionary representation of the amino acid properties."""
        return {
            "full_name": cls.full_name,
            "three_letter_code": cls.three_letter_code,
            "one_letter_code": cls.one_letter_code,
            "polar": cls.polar,
            "group": cls.group,
            "charge": cls.charge,
        }

    @classmethod
    def from_dict(cls, data):
        """Create an amino acid properties object from a dictionary."""
        return cls(
            full_name=data["full_name"],
            three_letter_code=data["three_letter_code"],
            one_letter_code=data["one_letter_code"],
            polar=data["polar"],
            group=data["group"],
            charge=data["charge"],
        )


class AminoAcids(Enum):
    """Class for managing amino acids."""

    # Defining aliases and their properties
    ALNANINE = ALA = A = AminoAcidProperties("alanine", "ALA", "A", False, "alkyl", "neutral")
    ARGININE = ARG = R = AminoAcidProperties("arginine", "ARG", "R", True, "basic", "positive")
    ASPARAGINE = ASN = N = AminoAcidProperties("asparagine", "ASN", "N", True, "neutral", "neutral")
    ASPARTIC_ACID = ASPARTATE = ASP = D = AminoAcidProperties(
        "aspartic acid", "ASP", "D", True, "acidic", "negative", "aspartate"
    )
    CYSTEINE = CYS = C = AminoAcidProperties("cysteine", "CYS", "C", True, "neutral", "neutral")
    GLUTAMINE = GLN = Q = AminoAcidProperties("glutamine", "GLN", "Q", True, "neutral", "neutral")
    GLUTAMIC_ACID = GLUTAMATE = GLU = E = AminoAcidProperties(
        "glutamic acid", "GLU", "E", True, "acidic", "negative", "glutamate"
    )
    GLYCINE = GLY = G = AminoAcidProperties("glycine", "GLY", "G", False, "neutral", "neutral")
    HISTIDINE = HIS = H = AminoAcidProperties("histidine", "HIS", "H", True, "basic", "positive")
    ISOLEUCINE = ILE = I = AminoAcidProperties("isoleucine", "ILE", "I", False, "alkyl", "neutral")  # noqa: E741
    LEUCINE = LEU = L = AminoAcidProperties("leucine", "LEU", "L", False, "alkyl", "neutral")
    LYSINE = LYS = K = AminoAcidProperties("lysine", "LYS", "K", True, "basic", "positive")
    METHIONINE = MET = M = AminoAcidProperties("methionine", "MET", "M", False, "neutral", "neutral")
    PHENYLALANINE = PHE = F = AminoAcidProperties("phenylalanine", "PHE", "F", False, "aromatic", "neutral")
    PROLINE = PRO = P = AminoAcidProperties("proline", "PRO", "P", False, "alkyl", "neutral")
    SERINE = SER = S = AminoAcidProperties("serine", "SER", "S", True, "neutral", "neutral")
    THREONINE = THR = T = AminoAcidProperties("threonine", "THR", "T", True, "neutral", "neutral")
    TRYPTOPHAN = TRP = W = AminoAcidProperties("tryptophan", "TRP", "W", True, "aromatic", "neutral")
    TYROSINE = TYR = Y = AminoAcidProperties("tyrosine", "TYR", "Y", False, "aromatic", "neutral")
    VALINE = VAL = V = AminoAcidProperties("valine", "VAL", "V", False, "alkyl", "neutral")

    def __repr__(self):
        """Return a string representation of the amino acid."""
        return f"{self.value.full_name} ({self.value.three_letter_code})"

    @property
    def full_name(self) -> str:
        """Return the full name of the amino acid."""
        return self.value.full_name

    @property
    def alt_name(self) -> str:
        """Return the alternative name of the amino acid."""
        if not self.value.alt_name:
            return self.value.full_name
        return self.value.alt_name

    @property
    def one_letter_code(self) -> str:
        """Return the one-letter code of the amino acid."""
        return self.value.one_letter_code

    @property
    def three_letter_code(self) -> str:
        """Return the three-letter code of the amino acid."""
        return self.value.three_letter_code

    @property
    def polar(self) -> bool:
        """Return whether the amino acid is polar."""
        return self.value.polar

    @property
    def group(self) -> str:
        """Return the group of the amino acid."""
        return self.value.group

    @property
    def charge(self) -> str:
        """Return the charge of the amino acid."""
        return self.value.charge

    @classmethod
    def filter_by_property(cls, **kwargs) -> list["AminoAcids"]:
        """Filter amino acids by given properties.

        Supported properties: polar (bool), group (str), charge (str).
        """
        return [aa for aa in cls if all(getattr(aa, prop) == value for prop, value in kwargs.items())]

    @classmethod
    def from_string(cls, value: str):
        """Convert a string to the corresponding AminoAcid enum member."""
        value = value.upper()

        # Check for matching names or aliases (one-letter, three-letter, full-name)
        for amino_acid in cls:
            if value.upper() in (
                amino_acid.name.upper(),
                amino_acid.value.one_letter_code.upper(),
                amino_acid.value.three_letter_code.upper(),
                amino_acid.value.full_name.upper(),
                amino_acid.value.alt_name.upper(),
            ):
                return amino_acid

        # If not found, raise an error
        msg = f"Invalid amino acid string: {value}"
        raise ValueError(msg)
