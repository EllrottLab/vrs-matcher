"""Unit tests for sample similarity metrics and ranking behavior.

These tests validate Jaccard computation, weighted concordance scoring,
pairwise matching, and one-vs-all ranking semantics.
"""

import pytest

from vrs_matcher.db import insert_alleles
from vrs_matcher.matcher import (
    MatchResult,
    jaccard,
    match_against_all,
    match_pair,
    weighted_concordance,
)
from vrs_matcher.models import GenotypeState, Zygosity

# S1: {aaa, bbb, ccc}  S2: {aaa, ddd}  S3: {aaa, bbb}
_ROWS = [
    ("S1", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 40.0, 30, None),
    ("S1", "ga4gh:VA.bbb", "1/1", "HOM_ALT", "chr2", 200, 50.0, 40, None),
    ("S1", "ga4gh:VA.ccc", "0/1", "HET", "chr3", 300, 45.0, 35, None),
    ("S2", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 35.0, 25, None),
    ("S2", "ga4gh:VA.ddd", "0/1", "HET", "chr4", 400, 38.0, 20, None),
    ("S3", "ga4gh:VA.aaa", "1/1", "HOM_ALT", "chr1", 100, 60.0, 50, None),
    ("S3", "ga4gh:VA.bbb", "0/1", "HET", "chr2", 200, 55.0, 45, None),
]


@pytest.fixture
def populated_db(db_conn):
    """Populate a temporary database with rows used by matcher tests.

    Args:
        db_conn: Open temporary SQLite connection fixture.

    Returns:
        sqlite3.Connection: The populated database connection.
    """

    insert_alleles(db_conn, _ROWS)
    return db_conn


class TestJaccard:
    """Tests for Jaccard set-overlap calculations."""

    def test_identical_sets(self):
        """Verify identical sets produce a Jaccard score of ``1.0``.

        Returns:
            None.
        """

        s = frozenset(["a", "b", "c"])
        assert jaccard(s, s) == pytest.approx(1.0)

    def test_disjoint_sets(self):
        """Verify disjoint sets produce a Jaccard score of ``0.0``.

        Returns:
            None.
        """

        assert jaccard(frozenset(["a"]), frozenset(["b"])) == pytest.approx(0.0)

    def test_partial_overlap(self):
        """Verify partial overlap produces the expected Jaccard ratio.

        Returns:
            None.
        """

        # |{a,b} ∩ {b,c}| / |{a,b,c}| = 1/3
        assert jaccard(frozenset(["a", "b"]), frozenset(["b", "c"])) == pytest.approx(1 / 3)

    def test_both_empty(self):
        """Verify both empty sets are treated as perfectly similar.

        Returns:
            None.
        """

        assert jaccard(frozenset(), frozenset()) == pytest.approx(1.0)

    def test_one_empty(self):
        """Verify one empty set yields zero similarity.

        Returns:
            None.
        """

        assert jaccard(frozenset(["a"]), frozenset()) == pytest.approx(0.0)


class TestWeightedConcordance:
    """Tests for per-variant concordance aggregation."""

    def test_same_genotype_string(self):
        """Verify identical genotype strings score full concordance.

        Returns:
            None.
        """

        het = GenotypeState(gt="0/1", zygosity=Zygosity.HET)
        assert weighted_concordance({"v1": het}, {"v1": het}) == pytest.approx(1.0)

    def test_same_vrs_id_different_zygosity(self):
        """Verify genotype disagreement on shared VRS IDs scores ``0.5``.

        Returns:
            None.
        """

        het = GenotypeState(gt="0/1", zygosity=Zygosity.HET)
        hom = GenotypeState(gt="1/1", zygosity=Zygosity.HOM_ALT)
        assert weighted_concordance({"v1": het}, {"v1": hom}) == pytest.approx(0.5)

    def test_no_call_scores_zero(self):
        """Verify NO_CALL in either sample yields zero concordance.

        Returns:
            None.
        """

        nc = GenotypeState(gt="./.", zygosity=Zygosity.NO_CALL)
        het = GenotypeState(gt="0/1", zygosity=Zygosity.HET)
        assert weighted_concordance({"v1": nc}, {"v1": het}) == pytest.approx(0.0)

    def test_no_shared_ids_returns_zero(self):
        """Verify no shared variants yields a concordance score of ``0.0``.

        Returns:
            None.
        """

        het = GenotypeState(gt="0/1", zygosity=Zygosity.HET)
        assert weighted_concordance({"v1": het}, {"v2": het}) == pytest.approx(0.0)

    def test_multiple_variants_averaged(self):
        """Verify concordance score is averaged across shared variants.

        Returns:
            None.
        """

        het = GenotypeState(gt="0/1", zygosity=Zygosity.HET)
        hom = GenotypeState(gt="1/1", zygosity=Zygosity.HOM_ALT)
        a = {"v1": het, "v2": het}
        b = {"v1": het, "v2": hom}
        # v1 → 1.0, v2 → 0.5 → avg = 0.75
        assert weighted_concordance(a, b) == pytest.approx(0.75)


class TestMatchPair:
    """Tests for pairwise sample matching results."""

    def test_overlap_metrics(self, populated_db):
        """Verify pairwise result includes expected overlap and totals.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        result = match_pair(populated_db, "S1", "S2")
        # S1={aaa,bbb,ccc} S2={aaa,ddd} → shared={aaa}, union=4
        assert result.jaccard == pytest.approx(1 / 4)
        assert "ga4gh:VA.aaa" in result.shared_vrs_ids
        assert result.total_a == 3
        assert result.total_b == 2

    def test_same_gt_gives_full_concordance(self, populated_db):
        """Verify shared equal genotypes produce full weighted concordance.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        # S1 and S2 share only aaa, both as "0/1" HET
        result = match_pair(populated_db, "S1", "S2")
        assert result.weighted_concordance == pytest.approx(1.0)

    def test_different_zygosity_gives_half_concordance(self, populated_db):
        """Verify differing genotypes on shared IDs yield half concordance.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        # S1(aaa=HET) vs S3(aaa=HOM_ALT): same VRS ID, different gt
        result = match_pair(populated_db, "S1", "S3")
        # shared = {aaa, bbb}; aaa→0.5, bbb→0.5 (S1 HOM_ALT vs S3 HET) → avg=0.5
        assert result.weighted_concordance == pytest.approx(0.5)

    def test_candidate_filter_raises_jaccard(self, populated_db):
        """Verify restricting to candidate IDs can raise similarity.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        candidates = frozenset(["ga4gh:VA.aaa"])
        result = match_pair(populated_db, "S1", "S2", candidate_vrs_ids=candidates)
        assert result.jaccard == pytest.approx(1.0)

    def test_result_fields(self, populated_db):
        """Verify match_pair returns a populated MatchResult object.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        result = match_pair(populated_db, "S1", "S2")
        assert isinstance(result, MatchResult)
        assert result.sample_a == "S1"
        assert result.sample_b == "S2"


class TestMatchAgainstAll:
    """Tests for one-vs-all sample ranking."""

    def test_sorted_by_descending_jaccard(self, populated_db):
        """Verify results are sorted by descending Jaccard score.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        results = match_against_all(populated_db, "S1")
        assert len(results) == 2
        assert results[0].jaccard >= results[1].jaccard

    def test_top_n_respected(self, populated_db):
        """Verify top-N truncation limits result count.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        results = match_against_all(populated_db, "S1", top_n=1)
        assert len(results) == 1

    def test_query_sample_not_compared_to_itself(self, populated_db):
        """Verify one-vs-all matching excludes the query sample itself.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        results = match_against_all(populated_db, "S1")
        for r in results:
            assert r.sample_b != "S1"

    def test_no_other_samples_returns_empty(self, db_conn):
        """Verify matching returns no results when only one sample exists.

        Args:
            db_conn: Open temporary SQLite connection fixture.

        Returns:
            None.
        """

        row = ("LONE", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 1, None, None, None)
        insert_alleles(db_conn, [row])
        assert match_against_all(db_conn, "LONE") == []

    def test_unknown_sample_raises_key_error(self, db_conn):
        """Verify match_against_all raises KeyError for an unregistered sample.

        Args:
            db_conn: Open temporary SQLite connection fixture.

        Returns:
            None.
        """

        import pytest

        with pytest.raises(KeyError, match="GHOST"):
            match_against_all(db_conn, "GHOST")


class TestMatchPairKeyError:
    """Tests that match_pair raises on unknown sample IDs."""

    def test_unknown_sample_a_raises(self, populated_db):
        """Verify KeyError is raised when sample_a is not in the index.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        import pytest

        with pytest.raises(KeyError, match="GHOST"):
            match_pair(populated_db, "GHOST", "S1")

    def test_unknown_sample_b_raises(self, populated_db):
        """Verify KeyError is raised when sample_b is not in the index.

        Args:
            populated_db: Populated temporary SQLite connection.

        Returns:
            None.
        """

        import pytest

        with pytest.raises(KeyError, match="GHOST"):
            match_pair(populated_db, "S1", "GHOST")
