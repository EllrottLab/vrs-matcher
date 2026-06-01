"""Plugin registration and discovery for matcher algorithms.

This module supports three plugin sources:

1. Built-in plugins registered by package code.
2. Third-party plugins exposed via Python entry points.
3. Local script plugins loaded from a file path.
"""

from __future__ import annotations

from importlib import metadata, util
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING, Any, Protocol

from .db import get_genotype_states, get_vrs_ids, list_samples, sample_exists

if TYPE_CHECKING:
    from .matcher import MatchResult
    from .models import GenotypeState

ENTRYPOINT_GROUP = "vrs_matcher.plugins"
PLUGIN_API_VERSION = "1"


class PluginError(RuntimeError):
    """Raised when plugin registration or resolution fails."""


class PluginContext:
    """Read-only database facade exposed to matcher plugins."""

    def __init__(self, conn) -> None:
        self.conn = conn

    def sample_exists(self, sample_id: str) -> bool:
        """Return whether ``sample_id`` is registered in the index."""

        return sample_exists(self.conn, sample_id)

    def get_vrs_ids(self, sample_id: str) -> frozenset[str]:
        """Return carried VRS IDs for ``sample_id``."""

        return get_vrs_ids(self.conn, sample_id)

    def get_genotype_states(self, sample_id: str) -> dict[str, GenotypeState]:
        """Return genotype-state mapping keyed by VRS ID for ``sample_id``."""

        return get_genotype_states(self.conn, sample_id)

    def list_samples(self) -> list[str]:
        """Return all registered sample IDs."""

        return list_samples(self.conn)


class MatcherPlugin(Protocol):
    """Structural contract for sample-matching plugins."""

    name: str
    api_version: str

    def match_pair(
        self,
        context: PluginContext,
        sample_a: str,
        sample_b: str,
        *,
        candidate_vrs_ids: frozenset[str] | None = None,
    ) -> MatchResult: ...

    def match_against_all(
        self,
        context: PluginContext,
        sample_id: str,
        *,
        top_n: int | None = None,
        candidate_vrs_ids: frozenset[str] | None = None,
    ) -> list[MatchResult]: ...


_BUILTIN_PLUGINS: dict[str, MatcherPlugin] = {}
_ENTRYPOINT_PLUGINS: dict[str, MatcherPlugin] | None = None


def _iter_entry_points() -> list[Any]:
    """Return entry points for the matcher plugin group."""

    all_eps = metadata.entry_points()
    if hasattr(all_eps, "select"):
        return list(all_eps.select(group=ENTRYPOINT_GROUP))
    return list(all_eps.get(ENTRYPOINT_GROUP, []))


def _instantiate_plugin(obj: Any) -> Any:
    """Instantiate plugin objects returned as classes or factories."""

    # Entry points may expose an instance, a plugin class (zero-arg ctor), or a
    # zero-arg factory returning either of those.
    if isinstance(obj, type):
        return obj()  # may raise TypeError; caller will handle

    if hasattr(obj, "match_pair") and hasattr(obj, "match_against_all"):
        return obj

    if callable(obj):
        return _instantiate_plugin(obj())

    return obj


def _validate_plugin(plugin: Any) -> MatcherPlugin:
    """Validate a plugin object and return it as ``MatcherPlugin``."""

    if not isinstance(getattr(plugin, "name", None), str) or not plugin.name:
        raise PluginError("Plugin must define a non-empty string 'name'.")

    api_version = getattr(plugin, "api_version", None)
    if api_version != PLUGIN_API_VERSION:
        raise PluginError(
            f"Plugin '{plugin.name}' API version '{api_version}' is not supported; "
            f"expected '{PLUGIN_API_VERSION}'."
        )

    if not callable(getattr(plugin, "match_pair", None)):
        raise PluginError(f"Plugin '{plugin.name}' must define match_pair(...).")
    if not callable(getattr(plugin, "match_against_all", None)):
        raise PluginError(f"Plugin '{plugin.name}' must define match_against_all(...).")

    return plugin


def register_builtin(plugin: MatcherPlugin, *, replace: bool = False) -> None:
    """Register a built-in plugin implementation."""

    validated = _validate_plugin(plugin)
    if validated.name in _BUILTIN_PLUGINS and not replace:
        raise PluginError(f"Built-in plugin already registered: {validated.name}")
    _BUILTIN_PLUGINS[validated.name] = validated


def clear_plugin_caches() -> None:
    """Reset cached plugin discovery state (useful in tests)."""

    global _ENTRYPOINT_PLUGINS
    _ENTRYPOINT_PLUGINS = None


def _load_entrypoint_plugins() -> dict[str, MatcherPlugin]:
    """Discover and cache entry-point plugins."""

    global _ENTRYPOINT_PLUGINS

    if _ENTRYPOINT_PLUGINS is not None:
        return _ENTRYPOINT_PLUGINS

    loaded: dict[str, MatcherPlugin] = {}
    for entry_point in _iter_entry_points():
        try:
            plugin_obj = _instantiate_plugin(entry_point.load())
            plugin = _validate_plugin(plugin_obj)
        except Exception:
            # Ignore broken third-party plugins; callers can still use built-ins.
            continue
        loaded.setdefault(plugin.name, plugin)

    _ENTRYPOINT_PLUGINS = loaded
    return loaded


def get_plugin(name: str) -> MatcherPlugin:
    """Resolve plugin by name from built-ins and entry points."""

    if name in _BUILTIN_PLUGINS:
        return _BUILTIN_PLUGINS[name]

    entrypoint_plugins = _load_entrypoint_plugins()
    if name in entrypoint_plugins:
        return entrypoint_plugins[name]

    raise KeyError(name)


def list_plugins() -> list[str]:
    """List available plugin names across all discovery sources."""

    names = set(_BUILTIN_PLUGINS)
    names.update(_load_entrypoint_plugins())
    return sorted(names)


def _load_module_from_path(path: Path) -> ModuleType:
    """Load a Python module from ``path``."""

    module_name = f"vrs_matcher_user_plugin_{abs(hash(path.resolve()))}"
    spec = util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise PluginError(f"Unable to load plugin module from: {path}")

    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_script_plugin(path: str | Path) -> MatcherPlugin:
    """Load a plugin from a local Python script.

    The script must expose ``create_plugin()`` returning a plugin instance.
    """

    module = _load_module_from_path(Path(path))
    factory = getattr(module, "create_plugin", None)
    if not callable(factory):
        raise PluginError("Plugin script must define create_plugin().")

    plugin = _instantiate_plugin(factory())
    return _validate_plugin(plugin)


def resolve_plugin(
    *,
    name: str | None = "identity",
    plugin_file: str | Path | None = None,
) -> MatcherPlugin:
    """Resolve a plugin by name or local script path."""

    if plugin_file is not None:
        plugin = load_script_plugin(plugin_file)
        # When using a script plugin, treat the default name ("identity") as
        # "no constraint" so callers don't have to pass name=None explicitly.
        if name not in (None, "identity") and plugin.name != name:
            raise KeyError(name)
        return plugin

    if name is None:
        raise PluginError("A plugin name is required when --plugin-file is not used.")

    return get_plugin(name)

    if name is None:
        raise PluginError("A plugin name is required when --plugin-file is not used.")

    return get_plugin(name)
