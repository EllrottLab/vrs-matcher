"""Example matcher plugin loaded from a local script.

This plugin computes the same identity metrics as the built-in algorithm but
applies a configurable floor to the reported Jaccard score.
"""

from __future__ import annotations

from vrs_matcher.matcher import MatchResult, jaccard, weighted_concordance
from vrs_matcher.plugins import PLUGIN_API_VERSION


class JaccardFloorPlugin:
    """Simple script plugin that enforces a minimum Jaccard score."""

    name = "jaccard-floor"
    api_version = PLUGIN_API_VERSION

    def __init__(self, floor: float = 0.25) -> None:
        self.floor = floor

    def match_pair(
        self,
        context,
        sample_a: str,
        sample_b: str,
        *,
        candidate_vrs_ids: frozenset[str] | None = None,
    ) -> MatchResult:
        for sid in (sample_a, sample_b):
            if not context.sample_exists(sid):
                raise KeyError(sid)

        ids_a = context.get_vrs_ids(sample_a)
        ids_b = context.get_vrs_ids(sample_b)
        states_a = context.get_genotype_states(sample_a)
        states_b = context.get_genotype_states(sample_b)

        if candidate_vrs_ids is not None:
            ids_a = ids_a & candidate_vrs_ids
            ids_b = ids_b & candidate_vrs_ids
            states_a = {k: v for k, v in states_a.items() if k in candidate_vrs_ids}
            states_b = {k: v for k, v in states_b.items() if k in candidate_vrs_ids}

        score = max(jaccard(ids_a, ids_b), self.floor)
        return MatchResult(
            sample_a=sample_a,
            sample_b=sample_b,
            jaccard=score,
            weighted_concordance=weighted_concordance(states_a, states_b),
            shared_vrs_ids=ids_a & ids_b,
            total_a=len(ids_a),
            total_b=len(ids_b),
        )

    def match_against_all(
        self,
        context,
        sample_id: str,
        *,
        top_n: int | None = None,
        candidate_vrs_ids: frozenset[str] | None = None,
    ) -> list[MatchResult]:
        if not context.sample_exists(sample_id):
            raise KeyError(sample_id)

        others = [sample for sample in context.list_samples() if sample != sample_id]
        results = [
            self.match_pair(
                context,
                sample_id,
                other,
                candidate_vrs_ids=candidate_vrs_ids,
            )
            for other in others
        ]
        results.sort(key=lambda result: result.jaccard, reverse=True)
        if top_n is not None:
            results = results[:top_n]
        return results


def create_plugin() -> JaccardFloorPlugin:
    """Factory required by vrs-matcher script plugin loading."""

    return JaccardFloorPlugin(floor=0.2)
