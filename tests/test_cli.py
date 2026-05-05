"""CLI tests for load/match/shared-variant command behavior.

These tests exercise command invocation, output formatting, and edge-case
handling for empty result scenarios.
"""

import pytest
from click.testing import CliRunner

from vrs_matcher.cli import cli
from vrs_matcher.db import insert_alleles, open_db

_ROWS = [
    ("S1", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 40.0, 30, None),
    ("S1", "ga4gh:VA.bbb", "1/1", "HOM_ALT", "chr2", 200, 50.0, 40, None),
    ("S2", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 35.0, 25, None),
]


@pytest.fixture
def test_db(tmp_path):
    """Create and populate a temporary database for CLI test execution.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        str: Filesystem path to the temporary SQLite database.
    """

    db = tmp_path / "test.db"
    conn = open_db(db)
    insert_alleles(conn, _ROWS)
    conn.close()
    return str(db)


def test_match_samples_output(test_db):
    """Verify match-samples prints expected summary fields and score.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    result = CliRunner().invoke(cli, ["match-samples", "S1", "S2", "--db", test_db])
    assert result.exit_code == 0
    assert "Jaccard" in result.output
    # S1={aaa,bbb} S2={aaa} → Jaccard=1/2=0.5
    assert "0.5000" in result.output
    assert "Shared variants:" in result.output


def test_match_samples_concordance(test_db):
    """Verify match-samples reports full concordance for matching genotypes.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    result = CliRunner().invoke(cli, ["match-samples", "S1", "S2", "--db", test_db])
    assert result.exit_code == 0
    # shared aaa: both 0/1 HET → concordance 1.0
    assert "1.0000" in result.output


def test_match_sample_against_all(test_db):
    """Verify match-sample returns at least one ranked peer sample.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    result = CliRunner().invoke(cli, ["match-sample", "S1", "--db", test_db, "--top", "5"])
    assert result.exit_code == 0
    assert "S2" in result.output


def test_match_sample_no_others(tmp_path):
    """Verify match-sample reports no peers when only one sample exists.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        None.
    """

    db = str(tmp_path / "lone.db")
    conn = open_db(db)
    insert_alleles(conn, [("LONE", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 1, None, None, None)])
    conn.close()
    result = CliRunner().invoke(cli, ["match-sample", "LONE", "--db", db])
    assert result.exit_code == 0
    assert "No other samples" in result.output


def test_shared_variants_lists_vrs_ids(test_db):
    """Verify shared-variants prints overlapping VRS identifiers.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    result = CliRunner().invoke(cli, ["shared-variants", "S1", "S2", "--db", test_db])
    assert result.exit_code == 0
    assert "ga4gh:VA.aaa" in result.output


def test_shared_variants_no_overlap(test_db):
    """Verify shared-variants prints nothing when overlap is empty.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    # S1 and a new sample with no shared variants
    conn = open_db(test_db)
    insert_alleles(conn, [("S3", "ga4gh:VA.zzz", "0/1", "HET", "chr9", 999, None, None, None)])
    conn.close()
    result = CliRunner().invoke(cli, ["shared-variants", "S2", "S3", "--db", test_db])
    assert result.exit_code == 0
    assert result.output.strip() == ""


def test_match_samples_unknown_sample_error(test_db):
    """Verify match-samples exits with an error for an unregistered sample.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    result = CliRunner().invoke(cli, ["match-samples", "GHOST", "S2", "--db", test_db])
    assert result.exit_code != 0
    assert "not found" in result.output.lower()


def test_match_sample_unknown_sample_error(test_db):
    """Verify match-sample exits with an error for an unregistered sample.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    result = CliRunner().invoke(cli, ["match-sample", "GHOST", "--db", test_db])
    assert result.exit_code != 0
    assert "not found" in result.output.lower()


def test_shared_variants_unknown_sample_error(test_db):
    """Verify shared-variants exits with an error for an unregistered sample.

    Args:
        test_db: Path to populated temporary SQLite database.

    Returns:
        None.
    """

    result = CliRunner().invoke(cli, ["shared-variants", "GHOST", "S2", "--db", test_db])
    assert result.exit_code != 0
    assert "not found" in result.output.lower()
