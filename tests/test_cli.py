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

# ── Identity-confirmation fixture rows ────────────────────────────────────────
# BASELINE (4 alleles) vs RESEQ (5 alleles, one zygosity mismatch) vs UNRELATED.
# Jaccard(BASELINE, RESEQ) = 4/5 = 0.8000
# Weighted concordance(BASELINE, RESEQ) = (1.0+1.0+1.0+0.5)/4 = 0.8750
_IDENTITY_ROWS = [
    ("BASELINE", "ga4gh:VA.v1", "0/1", "HET", "chr1", 101, 50.0, 30, "run1"),
    ("BASELINE", "ga4gh:VA.v2", "1/1", "HOM_ALT", "chr1", 102, 50.0, 30, "run1"),
    ("BASELINE", "ga4gh:VA.v3", "0/1", "HET", "chr1", 103, 50.0, 30, "run1"),
    ("BASELINE", "ga4gh:VA.v4", "1/1", "HOM_ALT", "chr1", 104, 50.0, 30, "run1"),
    ("RESEQ", "ga4gh:VA.v1", "0/1", "HET", "chr1", 101, 50.0, 30, "run2"),
    ("RESEQ", "ga4gh:VA.v2", "1/1", "HOM_ALT", "chr1", 102, 50.0, 30, "run2"),
    ("RESEQ", "ga4gh:VA.v3", "0/1", "HET", "chr1", 103, 50.0, 30, "run2"),
    ("RESEQ", "ga4gh:VA.v4", "0/1", "HET", "chr1", 104, 50.0, 30, "run2"),
    ("RESEQ", "ga4gh:VA.v5", "0/1", "HET", "chr1", 105, 50.0, 30, "run2"),
    ("UNRELATED", "ga4gh:VA.v1", "0/1", "HET", "chr1", 101, 50.0, 30, "run3"),
    ("UNRELATED", "ga4gh:VA.v6", "0/1", "HET", "chr1", 106, 50.0, 30, "run3"),
    ("UNRELATED", "ga4gh:VA.v7", "1/1", "HOM_ALT", "chr1", 107, 50.0, 30, "run3"),
]

# ── Release-QC fixture rows ───────────────────────────────────────────────────
# REFERENCE_RELEASE_SAMPLE (4 alleles) vs INCOMING_RELEASE_DUP (5 alleles,
# one zygosity mismatch) vs INCOMING_RELEASE_UNRELATED.
# Jaccard = 4/5 = 0.8000  |  Weighted concordance = 0.8750  (same pattern)
_RELEASE_QC_ROWS = [
    ("REFERENCE_RELEASE_SAMPLE", "ga4gh:VA.v1", "0/1", "HET", "chr1", 101, 50.0, 30, "reference_release"),
    ("REFERENCE_RELEASE_SAMPLE", "ga4gh:VA.v2", "1/1", "HOM_ALT", "chr1", 102, 50.0, 30, "reference_release"),
    ("REFERENCE_RELEASE_SAMPLE", "ga4gh:VA.v3", "0/1", "HET", "chr1", 103, 50.0, 30, "reference_release"),
    ("REFERENCE_RELEASE_SAMPLE", "ga4gh:VA.v4", "1/1", "HOM_ALT", "chr1", 104, 50.0, 30, "reference_release"),
    ("INCOMING_RELEASE_DUP", "ga4gh:VA.v1", "0/1", "HET", "chr1", 101, 50.0, 30, "incoming_release"),
    ("INCOMING_RELEASE_DUP", "ga4gh:VA.v2", "1/1", "HOM_ALT", "chr1", 102, 50.0, 30, "incoming_release"),
    ("INCOMING_RELEASE_DUP", "ga4gh:VA.v3", "0/1", "HET", "chr1", 103, 50.0, 30, "incoming_release"),
    ("INCOMING_RELEASE_DUP", "ga4gh:VA.v4", "0/1", "HET", "chr1", 104, 50.0, 30, "incoming_release"),
    ("INCOMING_RELEASE_DUP", "ga4gh:VA.v5", "0/1", "HET", "chr1", 105, 50.0, 30, "incoming_release"),
    ("INCOMING_RELEASE_UNRELATED", "ga4gh:VA.v1", "0/1", "HET", "chr1", 101, 50.0, 30, "incoming_release"),
    ("INCOMING_RELEASE_UNRELATED", "ga4gh:VA.v6", "0/1", "HET", "chr1", 106, 50.0, 30, "incoming_release"),
    ("INCOMING_RELEASE_UNRELATED", "ga4gh:VA.v7", "1/1", "HOM_ALT", "chr1", 107, 50.0, 30, "incoming_release"),
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


@pytest.fixture
def identity_db(tmp_path):
    """Populate a temporary database with identity-confirmation sample rows.

    Contains BASELINE, RESEQ, and UNRELATED samples.  Jaccard(BASELINE, RESEQ)
    = 0.8000; weighted concordance = 0.8750.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        str: Filesystem path to the temporary SQLite database.
    """

    db = tmp_path / "identity.db"
    conn = open_db(db)
    insert_alleles(conn, _IDENTITY_ROWS)
    conn.close()
    return str(db)


@pytest.fixture
def release_qc_db(tmp_path):
    """Populate a temporary database with release-QC sample rows.

    Contains REFERENCE_RELEASE_SAMPLE, INCOMING_RELEASE_DUP, and
    INCOMING_RELEASE_UNRELATED.  Jaccard(reference, dup) = 0.8000; weighted
    concordance = 0.8750.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        str: Filesystem path to the temporary SQLite database.
    """

    db = tmp_path / "release-qc.db"
    conn = open_db(db)
    insert_alleles(conn, _RELEASE_QC_ROWS)
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


def test_match_samples_identity_confirmation_scores(identity_db):
    """Verify CLI pairwise output for a near-duplicate identity-confirmation pair."""

    result = CliRunner().invoke(cli, ["match-samples", "BASELINE", "RESEQ", "--db", identity_db])
    assert result.exit_code == 0
    assert "Jaccard:" in result.output
    assert "0.8000" in result.output
    assert "Weighted concordance:" in result.output
    assert "0.8750" in result.output
    assert "Shared variants:" in result.output


def test_match_sample_identity_confirmation_top_hit(identity_db):
    """Verify one-vs-all ranking surfaces the resequenced sample as top hit."""

    result = CliRunner().invoke(
        cli, ["match-sample", "BASELINE", "--db", identity_db, "--top", "1"]
    )
    assert result.exit_code == 0

    data_lines = [line for line in result.output.splitlines() if line and not line.startswith("-")]
    assert len(data_lines) >= 2
    assert data_lines[1].startswith("RESEQ")


def test_match_samples_release_qc_scores(release_qc_db):
    """Verify CLI pairwise output for a release-style deduplication pair."""

    result = CliRunner().invoke(
        cli,
        ["match-samples", "REFERENCE_RELEASE_SAMPLE", "INCOMING_RELEASE_DUP", "--db", release_qc_db],
    )
    assert result.exit_code == 0
    assert "Jaccard:" in result.output
    assert "0.8000" in result.output
    assert "Weighted concordance:" in result.output
    assert "0.8750" in result.output
    assert "Shared variants:" in result.output


def test_match_sample_release_qc_top_hit(release_qc_db):
    """Verify the duplicate release sample ranks first in cohort screening."""

    result = CliRunner().invoke(
        cli, ["match-sample", "INCOMING_RELEASE_DUP", "--db", release_qc_db, "--top", "1"]
    )
    assert result.exit_code == 0

    data_lines = [line for line in result.output.splitlines() if line and not line.startswith("-")]
    assert len(data_lines) >= 2
    assert data_lines[1].startswith("REFERENCE_RELEASE_SAMPLE")
