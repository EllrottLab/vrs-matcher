"""Shared pytest fixtures for vrs_matcher test modules.

This module centralizes reusable fixtures so individual test files can focus on
behavior-specific assertions.
"""

import pytest

from vrs_matcher.db import open_db


@pytest.fixture
def db_conn(tmp_path):
    """Create a temporary SQLite connection with schema initialized.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Yields:
        sqlite3.Connection: Open database connection to a temporary test file.
    """

    conn = open_db(tmp_path / "test.db")
    yield conn
    conn.close()
