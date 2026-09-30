"""Data models for sample-level genotype projection.

The models in this module are intentionally minimal and focus on the fields
required for indexing and similarity scoring.
"""

from dataclasses import dataclass, field
from enum import StrEnum


class Zygosity(StrEnum):
    """Canonical zygosity categories used by the matcher.

    Attributes:
        HET: Heterozygous alternate call.
        HOM_ALT: Homozygous alternate call.
        REF: Homozygous reference call.
        NO_CALL: Missing or unknown genotype call.
    """

    HET = "HET"
    HOM_ALT = "HOM_ALT"
    REF = "REF"
    NO_CALL = "NO_CALL"


@dataclass
class GenotypeState:
    """Normalized genotype state for a sample at one VRS allele.

    Attributes:
        gt: Genotype string representation (for example, ``"0/1"`` or
            ``"1|1"``).
        zygosity: Normalized zygosity class inferred from GT alleles.
        depth: Optional read depth (DP) for this sample call.
        gq: Optional genotype quality (GQ) for this sample call.
    """

    gt: str
    zygosity: Zygosity
    depth: int | None = None
    gq: float | None = None


@dataclass
class SampleGenotype:
    """In-memory projection of all carried alleles for one sample.

    Attributes:
        sample_id: Unique sample identifier.
        alleles: Mapping of VRS allele identifier to genotype state.
    """

    sample_id: str
    alleles: dict[str, GenotypeState] = field(default_factory=dict)


@dataclass(frozen=True)
class KinshipResult:
    """KING equations 11 and 9 on jointly called markers; not identity scores."""

    sample_a: str
    sample_b: str
    panel_id: str
    n_common: int
    het_a: int
    het_b: int
    het_both: int
    opposite_hom: int
    kinship: float | None
    kinship_within_family: float | None
    ibs0_fraction: float | None
    status: str
    reason: str | None
    within_family_reason: str | None
    estimator: str = "KING-between-family-eq11"
    estimator_version: str = "1"


@dataclass
class KinshipMatches:
    """Scored top-N results and all unscorable peers, kept separate."""

    matches: list[KinshipResult]
    unscorable: list[KinshipResult]
