from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import pytest
from vrs_matcher.plugins import resolve_plugin, PluginContext
from vrs_matcher.db import open_db, insert_alleles


def _load_plugin(path: Path):
    spec = importlib.util.spec_from_file_location("_my_plugin_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load plugin: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "create_plugin"):
        raise RuntimeError("Plugin script must define create_plugin()")
    return module.create_plugin()


def test_my_plugin_differs_from_identity():
    """Test that my_plugin produces different scores than identity plugin.

    Creates synthetic data with varying allele prevalence to ensure
    weighted Jaccard differs from unweighted Jaccard.
    """
    plugin_file = Path("examples/plugins/my_plugin.py")
    if not plugin_file.exists():
        pytest.skip("examples/plugins/my_plugin.py not present in this environment")

    my_plugin = _load_plugin(plugin_file)

    # Create synthetic data with varying prevalence
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        db_path = tmp.name

    try:
        with open_db(db_path) as conn:
            # Create alleles with different prevalences:
            # - common_allele: present in 4/5 samples (80% prevalence, weight ~1.25)
            # - uncommon_allele: present in 2/5 samples (40% prevalence, weight = 2.5)
            # - rare_allele: present in 1/5 samples (20% prevalence, weight = 5.0)
            # This creates different weighted vs unweighted Jaccard scores
            rows = [
                # S1: has common + uncommon + rare
                ("S1", "ga4gh:VA.common", "0/1", "HET", "chr1", 100, 30.0, 20, None),
                ("S1", "ga4gh:VA.uncommon", "0/1", "HET", "chr1", 200, 30.0, 20, None),
                ("S1", "ga4gh:VA.rare", "0/1", "HET", "chr1", 300, 30.0, 20, None),

                # S2: has common + uncommon (overlaps with S1 on common+uncommon)
                ("S2", "ga4gh:VA.common", "0/1", "HET", "chr1", 100, 30.0, 20, None),
                ("S2", "ga4gh:VA.uncommon", "0/1", "HET", "chr1", 200, 30.0, 20, None),

                # S3: has common only
                ("S3", "ga4gh:VA.common", "0/1", "HET", "chr1", 100, 30.0, 20, None),

                # S4: has common only
                ("S4", "ga4gh:VA.common", "0/1", "HET", "chr1", 100, 30.0, 20, None),

                # S5: has common + unique allele (for diversity)
                ("S5", "ga4gh:VA.common", "0/1", "HET", "chr1", 100, 30.0, 20, None),
                ("S5", "ga4gh:VA.unique_to_s5", "0/1", "HET", "chr1", 400, 30.0, 20, None),
            ]
            insert_alleles(conn, rows)

        with open_db(db_path) as conn:
            ctx = PluginContext(conn)
            samples = ctx.list_samples()
            assert len(samples) >= 2, f"Expected at least 2 samples, got {len(samples)}"

            # Use S1 as the query sample (has common + uncommon + rare)
            sample_id = "S1"

            identity_plugin = resolve_plugin(name="identity")
            identity_results = identity_plugin.match_against_all(ctx, sample_id, top_n=10)
            custom_results = my_plugin.match_against_all(ctx, sample_id, top_n=10)

        by_identity = {r.sample_b: r for r in identity_results}
        by_custom = {r.sample_b: r for r in custom_results}
        common = set(by_identity) & set(by_custom)
        assert common, "No overlapping compared samples"

        # Check that at least one comparison differs
        differs = any(
            abs(float(by_custom[s].jaccard) - float(by_identity[s].jaccard)) > 1e-12
            for s in common
        )
        assert differs, (
            f"my-plugin scores are identical to identity; expected at least one difference. "
            f"Identity scores: {[(s, by_identity[s].jaccard) for s in common]}, "
            f"Custom scores: {[(s, by_custom[s].jaccard) for s in common]}"
        )

    finally:
        Path(db_path).unlink(missing_ok=True)
