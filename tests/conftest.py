"""Shared pytest fixtures for vrs_matcher test modules.

This module centralizes reusable fixtures so individual test files can focus on
behavior-specific assertions.
"""

import os

import pytest

from vrs_matcher.db import open_db


def pytest_addoption(parser: pytest.Parser) -> None:
    """Register opt-in CLI flag for integration tests."""

    parser.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="Run integration tests requiring network access and extra dependencies.",
    )


def pytest_configure(config: pytest.Config) -> None:
    """Declare custom markers used by this test suite."""

    config.addinivalue_line(
        "markers",
        "integration: slow tests requiring network access and optional dependencies",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip integration tests unless explicitly enabled."""

    allow_integration = config.getoption("--run-integration") or (
        os.getenv("RUN_INTEGRATION_TESTS") == "1"
    )
    if allow_integration:
        return

    skip_marker = pytest.mark.skip(
        reason=(
            "integration tests are disabled by default; pass --run-integration "
            "or set RUN_INTEGRATION_TESTS=1"
        )
    )
    for item in items:
        if "integration" in item.keywords:
            item.add_marker(skip_marker)


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
