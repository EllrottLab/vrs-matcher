"""Tests for matcher plugin registration and loading behavior."""

from __future__ import annotations

import inspect
from collections.abc import Iterator
from pathlib import Path
from typing import cast

import pytest

from vrs_matcher import plugins as plugin_mod
from vrs_matcher.db import insert_alleles
from vrs_matcher.matcher import match_against_all, match_pair
from vrs_matcher.plugins import (
    ENTRYPOINT_GROUP,
    PLUGIN_API_VERSION,
    MatcherPlugin,
    PluginError,
    clear_plugin_caches,
    list_plugins,
    resolve_plugin,
)

_ROWS = [
    ("S1", "ga4gh:VA.aaa", "0/1", "HET", "chr1", 100, 40.0, 30, None),
    ("S2", "ga4gh:VA.bbb", "0/1", "HET", "chr1", 110, 40.0, 30, None),
    ("S3", "ga4gh:VA.ccc", "0/1", "HET", "chr1", 120, 40.0, 30, None),
]


@pytest.fixture(autouse=True)
def reset_entrypoint_cache() -> Iterator[None]:
    """Ensure each test starts with a fresh entry-point plugin cache."""

    clear_plugin_caches()
    yield
    clear_plugin_caches()


def _write_script_plugin(path, *, name: str = "script-demo") -> str:
    path.write_text(
        f"""
from vrs_matcher.matcher import MatchResult
from vrs_matcher.plugins import PLUGIN_API_VERSION


class ScriptPlugin:
    name = \"{name}\"
    api_version = PLUGIN_API_VERSION

    def match_pair(self, context, sample_a, sample_b, *, candidate_vrs_ids=None):
        for sid in (sample_a, sample_b):
            if not context.sample_exists(sid):
                raise KeyError(sid)
        return MatchResult(
            sample_a=sample_a,
            sample_b=sample_b,
            jaccard=0.1234,
            weighted_concordance=0.0,
            shared_vrs_ids=frozenset(),
            total_a=0,
            total_b=0,
        )

    def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
        if not context.sample_exists(sample_id):
            raise KeyError(sample_id)
        others = [s for s in context.list_samples() if s != sample_id]
        rows = [
            MatchResult(
                sample_a=sample_id,
                sample_b=other,
                jaccard=0.0,
                weighted_concordance=0.0,
                shared_vrs_ids=frozenset(),
                total_a=0,
                total_b=0,
            )
            for other in sorted(others, reverse=True)
        ]
        return rows[:top_n] if top_n is not None else rows


def create_plugin():
    return ScriptPlugin()
""".strip(),
        encoding="utf-8",
    )
    return str(path)


def test_builtin_plugin_is_listed():
    """Built-in identity plugin should always be discoverable."""

    assert "identity" in list_plugins()


def test_clear_plugin_caches_resets_entrypoint_cache(monkeypatch):
    """The entry-point cache can be cleared and recomputed."""

    class FakeEntryPoint:
        def load(self):
            return object()

    class FakeEntryPoints(dict):
        def get(self, group, default=None):
            return [FakeEntryPoint()] if group == ENTRYPOINT_GROUP else default

    monkeypatch.setattr(plugin_mod.metadata, "entry_points", lambda: FakeEntryPoints())

    clear_plugin_caches()
    first = list_plugins()
    second = list_plugins()

    assert first == second
    assert "identity" in first


def test_entry_point_discovery_uses_legacy_get_branch(monkeypatch):
    """Legacy entry-point containers without select() should still work."""

    class LegacyPlugin:
        name = "legacy-demo"
        api_version = PLUGIN_API_VERSION

        def match_pair(self, context, sample_a, sample_b, *, candidate_vrs_ids=None):
            raise NotImplementedError

        def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
            raise NotImplementedError

    class FakeEntryPoint:
        def load(self):
            return LegacyPlugin

    class FakeEntryPoints(dict):
        def get(self, group, default=None):
            return [FakeEntryPoint()] if group == ENTRYPOINT_GROUP else default

    monkeypatch.setattr(plugin_mod.metadata, "entry_points", lambda: FakeEntryPoints())

    clear_plugin_caches()
    assert resolve_plugin(name="legacy-demo").name == "legacy-demo"


def test_broken_entry_point_is_ignored(monkeypatch):
    """Broken third-party entry points should not prevent built-ins from loading."""

    class GoodPlugin:
        name = "good-demo"
        api_version = PLUGIN_API_VERSION

        def match_pair(self, context, sample_a, sample_b, *, candidate_vrs_ids=None):
            raise NotImplementedError

        def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
            raise NotImplementedError

    def broken_factory():
        return object()

    class GoodEntryPoint:
        def load(self):
            return GoodPlugin

    class BrokenEntryPoint:
        def load(self):
            return broken_factory

    class FakeEntryPoints(dict):
        def select(self, *, group=None):
            return [BrokenEntryPoint(), GoodEntryPoint()] if group == ENTRYPOINT_GROUP else []

    monkeypatch.setattr(plugin_mod.metadata, "entry_points", lambda: FakeEntryPoints())

    clear_plugin_caches()
    names = list_plugins()
    assert "identity" in names
    assert "good-demo" in names
    assert "broken-demo" not in names


def test_script_plugin_can_drive_pairwise_matching(db_conn, tmp_path):
    """Script plugins can be loaded and used for match_pair execution."""

    insert_alleles(db_conn, _ROWS)
    plugin_file = _write_script_plugin(tmp_path / "my_plugin.py")

    result = match_pair(db_conn, "S1", "S2", plugin_file=plugin_file)
    assert result.jaccard == pytest.approx(0.1234)


def test_script_plugin_algorithm_name_mismatch_raises(db_conn, tmp_path):
    """Explicit algorithm names must match the plugin name for script plugins."""

    insert_alleles(db_conn, _ROWS)
    plugin_file = _write_script_plugin(tmp_path / "my_plugin.py", name="script-demo")

    with pytest.raises(PluginError, match="Matcher plugin not found"):
        match_pair(db_conn, "S1", "S2", algorithm="other", plugin_file=plugin_file)


def test_script_plugin_without_factory_raises(tmp_path):
    """Local plugin scripts must export create_plugin()."""

    plugin_file = tmp_path / "broken_plugin.py"
    plugin_file.write_text("NAME = 'broken'\n", encoding="utf-8")

    with pytest.raises(PluginError, match="create_plugin"):
        resolve_plugin(plugin_file=plugin_file)


def test_script_plugin_with_invalid_factory_raises(tmp_path):
    """A script factory returning an invalid object should fail validation."""

    plugin_file = tmp_path / "invalid_plugin.py"
    plugin_file.write_text(
        """
from vrs_matcher.plugins import PLUGIN_API_VERSION


def create_plugin():
    return type('InvalidPlugin', (), {'name': 'invalid', 'api_version': PLUGIN_API_VERSION})()
""".strip(),
        encoding="utf-8",
    )

    with pytest.raises(PluginError, match="match_pair"):
        resolve_plugin(plugin_file=plugin_file)


def test_script_plugin_controls_one_vs_all_order(db_conn, tmp_path):
    """Plugin-provided ranking order should be preserved by the wrapper."""

    insert_alleles(db_conn, _ROWS)
    plugin_file = _write_script_plugin(tmp_path / "my_plugin.py")

    results = match_against_all(db_conn, "S1", plugin_file=plugin_file)
    assert [result.sample_b for result in results] == ["S3", "S2"]


def test_example_plugin_in_examples_directory(db_conn):
    """The shipped example plugin should load and affect scores as documented."""

    insert_alleles(db_conn, _ROWS)
    plugin_file = Path("examples/plugins/jaccard_floor_plugin.py")

    result = match_pair(
        db_conn,
        "S1",
        "S2",
        algorithm="jaccard-floor",
        plugin_file=str(plugin_file),
    )

    assert plugin_file.exists()
    assert result.jaccard == pytest.approx(0.2)
    assert result.weighted_concordance == pytest.approx(0.0)

    ranked = match_against_all(
        db_conn,
        "S1",
        algorithm="jaccard-floor",
        plugin_file=str(plugin_file),
    )
    assert ranked
    assert all(match.jaccard >= 0.2 for match in ranked)


def test_entry_point_plugin_discovery(monkeypatch):
    """Entry-point plugins should be discoverable by name."""

    from vrs_matcher import plugins

    class EpPlugin:
        name = "ep-demo"
        api_version = PLUGIN_API_VERSION

        def match_pair(self, context, sample_a, sample_b, *, candidate_vrs_ids=None):
            raise NotImplementedError

        def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
            raise NotImplementedError

    class FakeEntryPoint:
        def load(self):
            return EpPlugin

    class FakeEntryPoints(list):
        def select(self, *, group=None):
            if group == ENTRYPOINT_GROUP:
                return self
            return []

    monkeypatch.setattr(
        plugins.metadata,
        "entry_points",
        lambda: FakeEntryPoints([FakeEntryPoint()]),
    )

    plugin = resolve_plugin(name="ep-demo")
    assert plugin.name == "ep-demo"


def test_register_builtin_rejects_duplicates(monkeypatch):
    """Built-in plugin names should be unique unless replacement is requested."""

    class TempPlugin:
        name = "temp-demo"
        api_version = PLUGIN_API_VERSION

        def match_pair(self, context, sample_a, sample_b, *, candidate_vrs_ids=None):
            raise NotImplementedError

        def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
            raise NotImplementedError

    snapshot = dict(plugin_mod._BUILTIN_PLUGINS)
    monkeypatch.setattr(plugin_mod, "_BUILTIN_PLUGINS", {})
    try:
        plugin_mod.register_builtin(cast(MatcherPlugin, TempPlugin()))
        with pytest.raises(PluginError, match="already registered"):
            plugin_mod.register_builtin(cast(MatcherPlugin, TempPlugin()))
    finally:
        monkeypatch.setattr(plugin_mod, "_BUILTIN_PLUGINS", snapshot)


def test_resolve_plugin_requires_name_without_plugin_file():
    """resolve_plugin should require a name if no script file is given."""

    with pytest.raises(PluginError, match="plugin name is required"):
        resolve_plugin(name=None)


def test_resolve_plugin_missing_name_raises_key_error():
    """Unknown plugin names should raise KeyError."""

    with pytest.raises(KeyError):
        resolve_plugin(name="definitely-missing")


def test_matcher_function_signatures_remain_notebook_friendly():
    """Notebook imports stay valid because core function shapes are preserved."""

    pair_params = list(inspect.signature(match_pair).parameters)
    against_all_params = list(inspect.signature(match_against_all).parameters)
    assert pair_params[:3] == ["conn", "sample_a", "sample_b"]
    assert against_all_params[:2] == ["conn", "sample_id"]
