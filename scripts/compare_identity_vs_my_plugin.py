#!/usr/bin/env python3
"""Compare IdentityMatcherPlugin vs my-plugin scores/ranks.

This script ensures the custom plugin is behaviorally different from identity
on a given database/query sample by checking that at least one compared pair
has a different primary score.

Usage:
  uv run python scripts/compare_identity_vs_my_plugin.py \
    --db matches.db \
    --sample SAMPLE_A \
    [--top 10]

Exit codes:
  0 -> differences observed
  1 -> no differences observed (or invalid input)
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

from vrs_matcher.storage import open_db

from vrs_matcher.matcher import match_sample


def _load_script_plugin(path: Path):
    spec = importlib.util.spec_from_file_location("_my_plugin_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load plugin script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "create_plugin"):
        raise RuntimeError("Plugin script must define create_plugin()")
    return module.create_plugin()


def _index_by_sample(results):
    return {r.sample_b: r for r in results}


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare identity vs my-plugin")
    parser.add_argument("--db", required=True, help="Path to vrs-matcher SQLite DB")
    parser.add_argument("--sample", required=True, help="Query sample ID")
    parser.add_argument("--top", type=int, default=10, help="Top-N matches to compare")
    parser.add_argument(
        "--plugin-file",
        default="examples/plugins/my_plugin.py",
        help="Path to my-plugin script file",
    )
    args = parser.parse_args()

    db_path = Path(args.db)
    plugin_path = Path(args.plugin_file)

    if not db_path.exists():
        print(f"ERROR: DB not found: {db_path}", file=sys.stderr)
        return 1
    if not plugin_path.exists():
        print(f"ERROR: Plugin script not found: {plugin_path}", file=sys.stderr)
        return 1

    my_plugin = _load_script_plugin(plugin_path)

    with open_db(str(db_path), read_only=True) as conn:
        rows = conn.execute(
            "SELECT sample_id FROM sample WHERE sample_id = ? LIMIT 1",
            (args.sample,),
        ).fetchall()
        if not rows:
            print(f"ERROR: sample not found: {args.sample}", file=sys.stderr)
            return 1

    identity_results = match_sample(
        str(db_path),
        args.sample,
        top_n=args.top,
        algorithm="identity",
    )

    my_results = match_sample(
        str(db_path),
        args.sample,
        top_n=args.top,
        plugin=my_plugin,
        algorithm=my_plugin.name,
    )

    idx_identity = _index_by_sample(identity_results)
    idx_my = _index_by_sample(my_results)
    common = sorted(set(idx_identity) & set(idx_my))

    if not common:
        print("ERROR: no overlapping matches to compare", file=sys.stderr)
        return 1

    print(f"Comparing sample: {args.sample}")
    print(f"DB: {db_path}")
    print(f"Top N: {args.top}")
    print(f"Plugin file: {plugin_path}")
    print()
    print("sample_b\tidentity_jaccard\tmy_score\tdelta")

    any_diff = False
    for sid in common:
        a = float(idx_identity[sid].jaccard)
        b = float(idx_my[sid].jaccard)
        d = b - a
        if abs(d) > 1e-12:
            any_diff = True
        print(f"{sid}\t{a:.6f}\t{b:.6f}\t{d:+.6f}")

    print()
    if any_diff:
        print("PASS: my-plugin differs from identity for at least one compared pair.")
        return 0

    print("FAIL: my-plugin did not differ from identity on compared pairs.", file=sys.stderr)
    print(
        "Hint: use a cohort where rare/common allele sharing differs between candidates.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
