"""Integration test validating population-structure signal from 1KGP data."""

from statistics import mean

import pytest

from vrs_matcher.db import open_db
from vrs_matcher.loader import load_samples
from vrs_matcher.matcher import match_against_all


@pytest.mark.integration
def test_intra_population_jaccard_exceeds_inter_population(
	tmp_path,
	annotated_vcf,
	sample_to_super_population,
) -> None:
	"""Assert same-super-population pairs are more similar than cross-population pairs."""

	db_path = tmp_path / "kgp.db"
	inserted = load_samples(annotated_vcf, db_path, source_dataset="1kg_chr22_16_17mb")
	assert inserted > 0

	conn = open_db(db_path)
	try:
		intra_scores: list[float] = []
		inter_scores: list[float] = []

		for sample_id, population in sample_to_super_population.items():
			for result in match_against_all(conn, sample_id):
				peer = result.sample_b
				peer_population = sample_to_super_population.get(peer)
				if peer_population is None:
					continue
				if peer_population == population:
					intra_scores.append(result.jaccard)
				else:
					inter_scores.append(result.jaccard)

		assert intra_scores, "No intra-population comparisons were produced"
		assert inter_scores, "No inter-population comparisons were produced"
		assert mean(intra_scores) > mean(inter_scores)
	finally:
		conn.close()

