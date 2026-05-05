"""Sample similarity scoring built on VRS allele identity.

This module provides pairwise and one-vs-all matching primitives based on VRS
set overlap and genotype concordance.
"""

from dataclasses import dataclass
from enum import StrEnum

from .db import get_genotype_states, get_vrs_ids, list_samples
from .models import GenotypeState, Zygosity


class MatchMode(StrEnum):
    """Named matching modes supported by the project model.

    Attributes:
        IDENTITY: Broad genome-wide similarity mode.
        RARE: Restriction to rare or high-impact variants.
        CANDIDATE: Restriction to user-provided candidate variants.
        PHENOTYPE: Restriction to phenotype-informed region/gene subsets.
    """

    IDENTITY = "identity"
    RARE = "rare"
    CANDIDATE = "candidate"
    PHENOTYPE = "phenotype"


@dataclass
class MatchResult:
    """Summary of sample-to-sample matching output.

    Attributes:
        sample_a: Query sample identifier.
        sample_b: Target sample identifier.
        jaccard: Jaccard similarity of carried VRS allele sets.
        weighted_concordance: Mean concordance score over shared VRS IDs.
        shared_vrs_ids: Set of VRS IDs carried by both samples.
        total_a: Number of carried VRS IDs in sample A.
        total_b: Number of carried VRS IDs in sample B.
    """

    sample_a: str
    sample_b: str
    jaccard: float
    weighted_concordance: float
    shared_vrs_ids: frozenset[str]
    total_a: int
    total_b: int


def jaccard(set_a: frozenset[str], set_b: frozenset[str]) -> float:
    """Compute Jaccard similarity between two VRS-ID sets.

    Jaccard is defined as ``|A ∩ B| / |A ∪ B|``.

    Args:
        set_a: First VRS-ID set.
        set_b: Second VRS-ID set.

    Returns:
        Jaccard coefficient in the range ``[0.0, 1.0]``.

    Notes:
        Returns ``1.0`` when both sets are empty.
    """
    if not set_a and not set_b:
        return 1.0
    union = len(set_a | set_b)
    return len(set_a & set_b) / union if union else 0.0


def _concordance(gt_a: GenotypeState, gt_b: GenotypeState) -> float:
    """Score genotype concordance for one shared VRS allele.

    Scoring rules:

    1. Same genotype string -> ``1.0``.
    2. Different genotype string -> ``0.5``.
    3. Either side is ``NO_CALL`` -> ``0.0``.

    Args:
        gt_a: Genotype state for sample A.
        gt_b: Genotype state for sample B.

    Returns:
        Concordance score for this shared allele.
    """
    if gt_a.zygosity == Zygosity.NO_CALL or gt_b.zygosity == Zygosity.NO_CALL:
        return 0.0
    if gt_a.gt == gt_b.gt:
        return 1.0
    return 0.5


def weighted_concordance(
    states_a: dict[str, GenotypeState],
    states_b: dict[str, GenotypeState],
) -> float:
    """Compute mean concordance across shared VRS IDs.

    Args:
        states_a: Mapping of VRS ID to genotype state for sample A.
        states_b: Mapping of VRS ID to genotype state for sample B.

    Returns:
        Mean concordance score across shared VRS IDs, or ``0.0`` when no
        variants are shared.
    """

    shared = set(states_a) & set(states_b)
    if not shared:
        return 0.0
    total = sum(_concordance(states_a[v], states_b[v]) for v in shared)
    return total / len(shared)


def match_pair(
    conn,
    sample_a: str,
    sample_b: str,
    *,
    candidate_vrs_ids: frozenset[str] | None = None,
) -> MatchResult:
    """Compute pairwise similarity between two samples.

    Args:
        conn: Open SQLite connection.
        sample_a: Query sample identifier.
        sample_b: Target sample identifier.
        candidate_vrs_ids: Optional restriction set of VRS IDs.

    Returns:
        A :class:`MatchResult` containing similarity metrics and summary counts.
    """

    ids_a = get_vrs_ids(conn, sample_a)
    ids_b = get_vrs_ids(conn, sample_b)

    if candidate_vrs_ids is not None:
        ids_a = ids_a & candidate_vrs_ids
        ids_b = ids_b & candidate_vrs_ids

    states_a = get_genotype_states(conn, sample_a)
    states_b = get_genotype_states(conn, sample_b)

    if candidate_vrs_ids is not None:
        states_a = {k: v for k, v in states_a.items() if k in candidate_vrs_ids}
        states_b = {k: v for k, v in states_b.items() if k in candidate_vrs_ids}

    return MatchResult(
        sample_a=sample_a,
        sample_b=sample_b,
        jaccard=jaccard(ids_a, ids_b),
        weighted_concordance=weighted_concordance(states_a, states_b),
        shared_vrs_ids=ids_a & ids_b,
        total_a=len(ids_a),
        total_b=len(ids_b),
    )


def match_against_all(
    conn,
    sample_id: str,
    *,
    top_n: int | None = None,
    candidate_vrs_ids: frozenset[str] | None = None,
) -> list[MatchResult]:
    """Match one sample against all other indexed samples.

    Args:
        conn: Open SQLite connection.
        sample_id: Query sample identifier.
        top_n: Optional maximum number of results to return.
        candidate_vrs_ids: Optional restriction set of VRS IDs.

    Returns:
        List of match results sorted by descending Jaccard score.
    """

    others = [s for s in list_samples(conn) if s != sample_id]
    results = [
        match_pair(conn, sample_id, other, candidate_vrs_ids=candidate_vrs_ids) for other in others
    ]
    results.sort(key=lambda r: r.jaccard, reverse=True)
    if top_n is not None:
        results = results[:top_n]
    return results
