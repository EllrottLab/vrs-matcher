"""Smoke tests for basic package import behavior."""

from vrs_matcher import __doc__


def test_package_importable() -> None:
    """Verify the package imports and exposes a module docstring.

    Returns:
        None.
    """

    assert __doc__ is not None
