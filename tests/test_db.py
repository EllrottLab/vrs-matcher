"""Unit tests for SQLite index helper functions.

These tests cover insert/query behavior for the sample-allele storage layer,
including upserts and optional field handling.
"""

import pytest

from vrs_matcher.db import (
    get_genotype_states,
    get_vrs_ids,
    insert_alleles,
    list_samples,
)
from vrs_matcher.models import Zygosity

_ROWS = [
    ("S1", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 40.0, 30, "cohort1"),
    ("S1", "ga4gh:VA.bbb", "1/1", "HOM_ALT", "chr2", 200, 50.0, 40, "cohort1"),
    ("S2", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 35.0, 25, "cohort1"),
    ("S2", "ga4gh:VA.ccc", "0/1", "HET", "chr3", 300, 45.0, 28, "cohort1"),
]


@pytest.fixture
def populated_db(db_conn):
    """Populate a temporary database with baseline sample-allele rows.

    Args:
        db_conn: Open temporary SQLite connection fixture.

    Returns:
        sqlite3.Connection: The populated database connection.
    """

    insert_alleles(db_conn, _ROWS)
    return db_conn


def test_insert_and_get_vrs_ids(populated_db):
    """Verify inserted VRS IDs are retrievable for a sample.

    Args:
        populated_db: Populated temporary SQLite connection.

    Returns:
        None.
    """

    assert get_vrs_ids(populated_db, "S1") == frozenset({"ga4gh:VA.aaa", "ga4gh:VA.bbb"})


def test_get_vrs_ids_unknown_sample_returns_empty(populated_db):
    """Verify unknown samples return an empty carried-VRS set.

    Args:
        populated_db: Populated temporary SQLite connection.

    Returns:
        None.
    """

    assert get_vrs_ids(populated_db, "UNKNOWN") == frozenset()


def test_get_genotype_states_gt(populated_db):
    """Verify genotype strings and zygosity values round-trip from storage.

    Args:
        populated_db: Populated temporary SQLite connection.

    Returns:
        None.
    """

    states = get_genotype_states(populated_db, "S1")
    assert states["ga4gh:VA.aaa"].gt == "0/1"
    assert states["ga4gh:VA.aaa"].zygosity == Zygosity.HET
    assert states["ga4gh:VA.bbb"].zygosity == Zygosity.HOM_ALT


def test_get_genotype_states_optional_fields(populated_db):
    """Verify optional QC fields are retrieved with expected values.

    Args:
        populated_db: Populated temporary SQLite connection.

    Returns:
        None.
    """

    states = get_genotype_states(populated_db, "S1")
    assert states["ga4gh:VA.aaa"].gq == pytest.approx(40.0)
    assert states["ga4gh:VA.aaa"].depth == 30


def test_list_samples(populated_db):
    """Verify sample listing returns distinct IDs in sorted order.

    Args:
        populated_db: Populated temporary SQLite connection.

    Returns:
        None.
    """

    assert list_samples(populated_db) == ["S1", "S2"]


def test_upsert_replaces_existing_row(db_conn):
    """Verify insert uses replace semantics for duplicate primary keys.

    Args:
        db_conn: Open temporary SQLite connection fixture.

    Returns:
        None.
    """

    insert_alleles(db_conn, [("S1", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 40.0, 30, None)])
    insert_alleles(db_conn, [("S1", "ga4gh:VA.aaa", "1/1", "HOM_ALT", "chr1", 100, 50.0, 35, None)])
    states = get_genotype_states(db_conn, "S1")
    assert states["ga4gh:VA.aaa"].zygosity == Zygosity.HOM_ALT


def test_null_optional_fields(db_conn):
    """Verify nullable QC columns are mapped back to ``None`` in model state.

    Args:
        db_conn: Open temporary SQLite connection fixture.

    Returns:
        None.
    """

    insert_alleles(db_conn, [("S1", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, None, None, None)])
    states = get_genotype_states(db_conn, "S1")
    assert states["ga4gh:VA.aaa"].gq is None
    assert states["ga4gh:VA.aaa"].depth is None
