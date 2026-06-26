from __future__ import annotations

from collections.abc import Mapping, Set

from vrs_matcher.matcher import MatchResult, weighted_concordance
from vrs_matcher.plugins import PLUGIN_API_VERSION


class MyPlugin:
    """Rare-variant-weighted matcher.

    Differences from IdentityMatcherPlugin:
    - Uses a weighted overlap score (stored in MatchResult.jaccard) where rarer
      alleles contribute more than common alleles.
    - Uses cohort prevalence (fraction of samples containing each allele) as a
      rarity proxy derived from ``context.list_samples()`` + ``context.get_vrs_ids``.
    - Keeps weighted_concordance unchanged for compatibility.
    """

    name = "my-plugin"
    api_version = PLUGIN_API_VERSION

    # Simple, explainable weighting with caps for stability.
    # prevalence in (0, 1] -> weight = 1 / prevalence, then clipped
    min_weight = 1.0
    max_weight = 10.0

    def _build_weight_map(self, context, candidate_vrs_ids: Set[str] | None) -> dict[str, float]:
        """Compute per-allele rarity weights from cohort prevalence.

        prevalence(vrs_id) = count(samples containing vrs_id) / total_samples
        weight(vrs_id) = clip(1 / prevalence, min_weight, max_weight)
        """
        samples = list(context.list_samples())
        if not samples:
            return {}

        counts: dict[str, int] = {}
        for sid in samples:
            ids = set(context.get_vrs_ids(sid))
            if candidate_vrs_ids is not None:
                ids &= set(candidate_vrs_ids)
            for vid in ids:
                counts[vid] = counts.get(vid, 0) + 1

        n = float(len(samples))
        weights: dict[str, float] = {}
        for vid, c in counts.items():
            prevalence = c / n
            # prevalence should never be 0 here, but guard defensively.
            raw = 1.0 / prevalence if prevalence > 0 else self.max_weight
            if raw < self.min_weight:
                raw = self.min_weight
            elif raw > self.max_weight:
                raw = self.max_weight
            weights[vid] = raw
        return weights

    @staticmethod
    def _weighted_jaccard(ids_a: Set[str], ids_b: Set[str], weights: Mapping[str, float]) -> float:
        """Weighted Jaccard on sets using per-id weight map.

        score = sum(w_i for i in A∩B) / sum(w_i for i in A∪B)
        Defaults missing weights to 1.0.
        """
        union = ids_a | ids_b
        if not union:
            return 0.0

        inter = ids_a & ids_b
        num = sum(float(weights.get(v, 1.0)) for v in inter)
        den = sum(float(weights.get(v, 1.0)) for v in union)
        return num / den if den else 0.0

    def match_pair(self, context, sample_a, sample_b, *, candidate_vrs_ids=None):
        for sid in (sample_a, sample_b):
            if not context.sample_exists(sid):
                raise KeyError(sid)

        ids_a = set(context.get_vrs_ids(sample_a))
        ids_b = set(context.get_vrs_ids(sample_b))
        states_a = dict(context.get_genotype_states(sample_a))
        states_b = dict(context.get_genotype_states(sample_b))

        if candidate_vrs_ids is not None:
            candidate_vrs_ids = set(candidate_vrs_ids)
            ids_a &= candidate_vrs_ids
            ids_b &= candidate_vrs_ids
            states_a = {k: v for k, v in states_a.items() if k in candidate_vrs_ids}
            states_b = {k: v for k, v in states_b.items() if k in candidate_vrs_ids}

        weights = self._build_weight_map(context, candidate_vrs_ids)
        weighted_overlap = self._weighted_jaccard(ids_a, ids_b, weights)

        return MatchResult(
            sample_a=sample_a,
            sample_b=sample_b,
            # Stored in jaccard field for CLI compatibility.
            jaccard=weighted_overlap,
            weighted_concordance=weighted_concordance(states_a, states_b),
            shared_vrs_ids=ids_a & ids_b,
            total_a=len(ids_a),
            total_b=len(ids_b),
        )

    def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
        if not context.sample_exists(sample_id):
            raise KeyError(sample_id)

        results = [
            self.match_pair(context, sample_id, other, candidate_vrs_ids=candidate_vrs_ids)
            for other in context.list_samples()
            if other != sample_id
        ]
        results.sort(key=lambda result: result.jaccard, reverse=True)
        return results[:top_n] if top_n is not None else results


def create_plugin():
    return MyPlugin()
