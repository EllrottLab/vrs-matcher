"""SQLite storage helpers for the sample-allele index.

This module owns schema creation and provides small query/update helpers used by
the loader, matcher, and CLI layers.
"""

import sqlite3
from collections.abc import Iterable
from pathlib import Path

from .models import GenotypeState, Zygosity

_SCHEMA = """
CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT PRIMARY KEY
);
CREATE TABLE IF NOT EXISTS sample_allele (
    sample_id      TEXT,
    vrs_id         TEXT,
    gt             TEXT,
    zygosity       TEXT,
    chrom          TEXT,
    pos            INTEGER,
    gq             REAL,
    dp             INTEGER,
    source_dataset TEXT,
    PRIMARY KEY (sample_id, vrs_id)
);
CREATE INDEX IF NOT EXISTS idx_vrs_id    ON sample_allele(vrs_id);
CREATE INDEX IF NOT EXISTS idx_sample_id ON sample_allele(sample_id);
"""


def open_db(path: str | Path) -> sqlite3.Connection:
    """Open or create the sample-allele index database.

    Args:
        path: Filesystem path to the SQLite database file.

    Returns:
        An open SQLite connection with row factory configured and schema
        initialized.
    """

    conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    conn.executescript(_SCHEMA)
    conn.commit()
    return conn


def register_samples(conn: sqlite3.Connection, sample_ids: Iterable[str]) -> None:
    """Register sample IDs in the samples table.

    Idempotent: already-present IDs are silently ignored.  Call this with the
    full VCF-header sample list during ingestion so that samples with zero
    surviving alleles are still discoverable.

    Args:
        conn: Open SQLite connection.
        sample_ids: Iterable of sample identifiers to register.

    Returns:
        None.
    """

    conn.executemany(
        "INSERT OR IGNORE INTO samples (sample_id) VALUES (?)",
        ((sid,) for sid in sample_ids),
    )
    conn.commit()


def sample_exists(conn: sqlite3.Connection, sample_id: str) -> bool:
    """Return whether a sample ID is present in the index.

    A sample is considered present if it was registered during ingestion, even
    if it has zero surviving allele rows after filtering.

    Args:
        conn: Open SQLite connection.
        sample_id: Sample identifier to look up.

    Returns:
        ``True`` if the sample is registered, ``False`` otherwise.
    """

    cur = conn.execute("SELECT 1 FROM samples WHERE sample_id = ?", (sample_id,))
    return cur.fetchone() is not None


def insert_alleles(conn: sqlite3.Connection, rows: Iterable[tuple]) -> None:
    """Insert or replace allele rows into the index.

    Also registers the sample IDs found in ``rows`` so that callers which
    bypass the VCF loader (e.g. tests) are still reflected in
    :func:`list_samples` and :func:`sample_exists`.

    Args:
        conn: Open SQLite connection.
        rows: Iterable of row tuples matching the ``sample_allele`` column
            order:
            ``(sample_id, vrs_id, gt, zygosity, chrom, pos, gq, dp, source_dataset)``.

    Returns:
        None.
    """

    rows_list = list(rows)
    register_samples(conn, dict.fromkeys(row[0] for row in rows_list))
    conn.executemany(
        """
        INSERT OR REPLACE INTO sample_allele
            (sample_id, vrs_id, gt, zygosity, chrom, pos, gq, dp, source_dataset)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows_list,
    )
    conn.commit()


def get_vrs_ids(conn: sqlite3.Connection, sample_id: str) -> frozenset[str]:
    """Return the set of VRS IDs carried by a sample.

    Args:
        conn: Open SQLite connection.
        sample_id: Sample identifier to query.

    Returns:
        A frozenset of VRS IDs associated with ``sample_id``.
    """

    cur = conn.execute("SELECT vrs_id FROM sample_allele WHERE sample_id = ?", (sample_id,))
    return frozenset(row["vrs_id"] for row in cur)


def get_genotype_states(conn: sqlite3.Connection, sample_id: str) -> dict[str, GenotypeState]:
    """Return genotype states keyed by VRS ID for one sample.

    Args:
        conn: Open SQLite connection.
        sample_id: Sample identifier to query.

    Returns:
        A dictionary mapping each VRS ID to its corresponding
        :class:`vrs_matcher.models.GenotypeState`.
    """

    cur = conn.execute(
        "SELECT vrs_id, gt, zygosity, gq, dp FROM sample_allele WHERE sample_id = ?",
        (sample_id,),
    )
    return {
        row["vrs_id"]: GenotypeState(
            gt=row["gt"],
            zygosity=Zygosity(row["zygosity"]),
            gq=row["gq"],
            depth=row["dp"],
        )
        for row in cur
    }


def list_samples(conn: sqlite3.Connection) -> list[str]:
    """Return all sample IDs registered in the index.

    Includes samples that were processed during ingestion but had zero
    surviving allele rows after filtering.

    Args:
        conn: Open SQLite connection.

    Returns:
        Sorted list of unique sample identifiers.
    """

    cur = conn.execute("SELECT sample_id FROM samples ORDER BY sample_id")
    return [row["sample_id"] for row in cur]
