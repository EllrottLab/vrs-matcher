"""Manichaikul et al. (2010), equations 11 and 9 on pairwise complete calls."""

import heapq
import itertools

from .models import KinshipMatches, KinshipResult
from .plugins import PLUGIN_API_VERSION, PluginContext, PluginError

_GENOTYPE_BATCH_SIZE = 128


class _ReverseLexical:
    """Reverse lexical ordering so the heap root is the worst tied sample ID."""

    def __init__(self, value: str) -> None:
        self.value = value

    def __lt__(self, other: "_ReverseLexical") -> bool:
        return self.value > other.value


def _score(sample_a, sample_b, panel_id, calls_a, calls_b, candidates):
    common = calls_a.keys() & calls_b.keys()
    if candidates is not None:
        common &= candidates
    m = len(common)
    ha = hb = hh = opposite = 0
    for marker in common:
        a, b = calls_a[marker], calls_b[marker]
        ha += a == 1
        hb += b == 1
        hh += a == b == 1
        opposite += abs(a - b) == 2
    reason = "no_common_calls" if not m else "zero_heterozygotes" if min(ha, hb) == 0 else None
    return KinshipResult(
        sample_a=sample_a,
        sample_b=sample_b,
        panel_id=panel_id,
        n_common=m,
        het_a=ha,
        het_b=hb,
        het_both=hh,
        opposite_hom=opposite,
        kinship=0.5 - (ha + hb - 2 * hh + 4 * opposite) / (4 * min(ha, hb))
        if min(ha, hb)
        else None,
        kinship_within_family=(hh - 2 * opposite) / (ha + hb) if ha + hb else None,
        ibs0_fraction=opposite / m if m else None,
        status="unscorable" if reason else "ok",
        reason=reason,
        within_family_reason=None if ha + hb else reason,
    )


class KingRobustPlugin:
    name = "king-robust"
    api_version = PLUGIN_API_VERSION

    def _panel(self, context, candidates):
        panel = context.get_genotype_panel()
        if candidates is not None and candidates - context.get_panel_ids():
            raise PluginError("Candidate VRS IDs must all belong to the KING panel.")
        return panel["panel_id"]

    def match_pair(self, context: PluginContext, sample_a, sample_b, *, candidate_vrs_ids=None):
        panel_id = self._panel(context, candidate_vrs_ids)
        return _score(
            sample_a,
            sample_b,
            panel_id,
            context.get_called_genotypes(sample_a),
            context.get_called_genotypes(sample_b),
            candidate_vrs_ids,
        )

    def iter_all_pairs(
        self, context: PluginContext, *, candidate_vrs_ids=None, batch_size=_GENOTYPE_BATCH_SIZE
    ):
        """Yield distinct pair scores while retaining at most two call batches."""
        if batch_size < 1:
            raise PluginError("batch_size must be positive.")
        panel_id = self._panel(context, candidate_vrs_ids)
        samples = context.list_samples()
        for left_start in range(0, len(samples), batch_size):
            left_ids = samples[left_start : left_start + batch_size]
            calls_left = context.get_called_genotypes_many(left_ids, batch_size=batch_size)
            for right_start in range(left_start, len(samples), batch_size):
                if right_start == left_start:
                    if len(left_ids) < 2:
                        continue
                    for sample_a, sample_b in itertools.combinations(left_ids, 2):
                        yield _score(
                            sample_a,
                            sample_b,
                            panel_id,
                            calls_left[sample_a],
                            calls_left[sample_b],
                            candidate_vrs_ids,
                        )
                    continue
                right_ids = samples[right_start : right_start + batch_size]
                calls_right = context.get_called_genotypes_many(right_ids, batch_size=batch_size)
                for sample_a in left_ids:
                    for sample_b in right_ids:
                        yield _score(
                            sample_a,
                            sample_b,
                            panel_id,
                            calls_left[sample_a],
                            calls_right[sample_b],
                            candidate_vrs_ids,
                        )
                del calls_right

    def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
        if top_n is not None and top_n < 1:
            raise PluginError("top_n must be positive.")
        panel_id = self._panel(context, candidate_vrs_ids)
        samples = context.list_samples()
        if sample_id not in samples:
            raise KeyError(sample_id)
        peer_ids = [peer for peer in samples if peer != sample_id]
        query = context.get_called_genotypes(sample_id)
        matches, unscorable = [], []
        heap = []
        for start in range(0, len(peer_ids), _GENOTYPE_BATCH_SIZE):
            peer_batch = peer_ids[start : start + _GENOTYPE_BATCH_SIZE]
            calls = context.get_called_genotypes_many(peer_batch)
            for peer in peer_batch:
                result = _score(
                    sample_id,
                    peer,
                    panel_id,
                    query,
                    calls[peer],
                    candidate_vrs_ids,
                )
                if result.kinship is None:
                    unscorable.append(result)
                elif top_n is None:
                    matches.append(result)
                else:
                    entry = (result.kinship, _ReverseLexical(result.sample_b), result)
                    if len(heap) < top_n:
                        heapq.heappush(heap, entry)
                    elif heap[0][:2] < entry[:2]:
                        heapq.heapreplace(heap, entry)
        if top_n is not None:
            matches = [entry[2] for entry in heap]
        matches.sort(key=lambda r: (-r.kinship, r.sample_b))
        return KinshipMatches(matches=matches[:top_n], unscorable=unscorable)
