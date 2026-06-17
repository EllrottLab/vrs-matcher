from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from vrs_matcher.matcher import match_sample
from vrs_matcher.storage import open_db


def _load_plugin(path: Path):
    spec = importlib.util.spec_from_file_location("_my_plugin_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load plugin: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "create_plugin"):
        raise RuntimeError("Plugin script must define create_plugin()")
    return module.create_plugin()


@pytest.mark.integration
def test_my_plugin_differs_from_identity():
    db = Path("tests/data/matches.db")
    plugin_file = Path("examples/plugins/my_plugin.py")

    if not db.exists():
        pytest.skip("tests/data/matches.db not present in this environment")
    if not plugin_file.exists():
        pytest.skip("examples/plugins/my_plugin.py not present in this environment")

    plugin = _load_plugin(plugin_file)

    # Pick a query sample by reading one from DB.
    # Uses only public APIs available in this repo.
    with open_db(str(db), read_only=True) as conn:
        rows = conn.execute("SELECT sample_id FROM sample ORDER BY sample_id LIMIT 1").fetchall()
        if not rows:
            pytest.skip("No samples found in tests/data/matches.db")
        query = rows[0][0]

    identity = match_sample(str(db), query, top_n=10, algorithm="identity")
    custom = match_sample(str(db), query, top_n=10, plugin=plugin, algorithm=plugin.name)

    by_identity = {r.sample_b: r for r in identity}
    by_custom = {r.sample_b: r for r in custom}
    common = set(by_identity) & set(by_custom)

    assert common, "No overlapping compared samples"

    differs = any(
        abs(float(by_custom[s].jaccard) - float(by_identity[s].jaccard)) > 1e-12
        for s in common
    )
    assert differs, "my-plugin scores are identical to identity; expected at least one difference"
