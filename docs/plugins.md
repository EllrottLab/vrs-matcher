# Matcher plugins

`vrs-matcher` can run matching algorithms supplied by users, not just the
built-in `identity` algorithm. Write a plugin when the built-in matcher is
close, but not quite the scoring or ranking logic you need — for example:

- weighting rare or clinically important variants more heavily,
- restricting matching to a candidate gene list or panel,
- down-weighting noisy loci for cohort deduplication, or
- penalizing discordant zygosity at high-confidence loci more strongly.

If your only goal is to restrict the comparison to a subset of variants, first
check whether the `candidate_vrs_ids` argument is enough. A plugin is the right
tool when you need different scoring or ranking behavior.

## What a plugin can and cannot change

A plugin changes **how samples are scored and ranked** after data are loaded
into the SQLite index. It operates on indexed sample-level allele/genotype
data, not raw VCF records.

Plugins do **not** replace VCF parsing, VRS allele extraction, SQLite storage,
or the `MatchResult` output structure.

## Plugin sources and selection

There are three plugin sources, selected on the command line:

| Source | How to select | Notes |
| --- | --- | --- |
| **Built-in** (shipped with `vrs-matcher`) | `--algorithm identity` | |
| **Entry-point** (from an installed package) | `--algorithm my-plugin` | See [Packaging](#packaging-a-reusable-plugin) |
| **Local script** (a `.py` file) | `--plugin-file my_plugin.py` | Usually no `--algorithm` needed |

When using `--plugin-file`, you do not need `--algorithm`; if you supply it
anyway, it must match the plugin's `name`.

The recommended path is to start with a local script plugin, test it on your
cohort, then package it as an entry-point plugin if you want to share it.

### Loading order

At runtime, `resolve_plugin(...)` (called by the public functions in
`src/vrs_matcher/matcher.py`) selects a plugin in this order:

1. **Built-in plugins** — registered with `register_builtin(...)` when
   `vrs_matcher.matcher` is imported (the identity matcher lives there).
2. **Package plugins** — discovered from the `vrs_matcher.plugins` entry-point
   group via package metadata.
3. **Script plugins** — loaded from the `--plugin-file` path; the file must
   define `create_plugin()`.

## Quick start

```bash
# See available built-in and installed plugins
uv run vrs-matcher plugins list

# Run the built-in algorithm explicitly
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --algorithm identity

# Run a local script plugin
uv run vrs-matcher match-sample SAMPLE_A --db matches.db \
  --plugin-file examples/plugins/jaccard_floor_plugin.py
```

`--plugin-file` also works with `match-samples` (pairwise) and
`shared-variants`:

```bash
uv run vrs-matcher match-samples SAMPLE_A SAMPLE_B --db matches.db \
  --plugin-file examples/plugins/jaccard_floor_plugin.py
```

The repository includes a complete worked example,
`examples/plugins/jaccard_floor_plugin.py`. It behaves like the built-in
identity matcher but enforces a minimum Jaccard score floor — a teaching
example, not a recommended production method.

## Plugin contract

A plugin object must define:

- `name` — unique plugin name, e.g. `"identity"` or `"rare-priority"`
- `api_version` — currently `"1"` (use `PLUGIN_API_VERSION`)
- `match_pair(context, sample_a, sample_b, *, candidate_vrs_ids=None) -> MatchResult`
- `match_against_all(context, sample_id, *, top_n=None, candidate_vrs_ids=None) -> list[MatchResult]`

A **script plugin** file must additionally define `create_plugin()`, returning
the plugin object.

### Input: the `context` object

Plugins receive a read-only `context` object instead of direct SQL access, so
plugin code stays decoupled from the database schema. By the time a plugin
runs, samples are loaded, VRS IDs are normalized, and genotype states are
available per sample — the plugin only needs to focus on scoring.

Available helpers:

- `context.sample_exists(sample_id)`
- `context.get_vrs_ids(sample_id)`
- `context.get_genotype_states(sample_id)`
- `context.list_samples()`

### Output: the `MatchResult`

Plugins return `MatchResult` objects (from `vrs_matcher.matcher`), populating:
`sample_a`, `sample_b`, `jaccard`, `weighted_concordance`, `shared_vrs_ids`,
`total_a`, `total_b`.

> **Note:** the CLI and notebook interfaces report the primary score from the
> `jaccard` field. If your algorithm is not truly Jaccard-based and you reuse
> that field for a custom score, document it clearly for users.

## Minimal script plugin template

Copy this into `my_plugin.py` and modify the scoring logic.

```python
from vrs_matcher.matcher import MatchResult, jaccard, weighted_concordance
from vrs_matcher.plugins import PLUGIN_API_VERSION


class MyPlugin:
    name = "my-plugin"
    api_version = PLUGIN_API_VERSION

    def match_pair(self, context, sample_a, sample_b, *, candidate_vrs_ids=None):
        for sid in (sample_a, sample_b):
            if not context.sample_exists(sid):
                raise KeyError(sid)

        ids_a = context.get_vrs_ids(sample_a)
        ids_b = context.get_vrs_ids(sample_b)
        states_a = context.get_genotype_states(sample_a)
        states_b = context.get_genotype_states(sample_b)

        if candidate_vrs_ids is not None:
            ids_a = ids_a & candidate_vrs_ids
            ids_b = ids_b & candidate_vrs_ids
            states_a = {k: v for k, v in states_a.items() if k in candidate_vrs_ids}
            states_b = {k: v for k, v in states_b.items() if k in candidate_vrs_ids}

        return MatchResult(
            sample_a=sample_a,
            sample_b=sample_b,
            jaccard=jaccard(ids_a, ids_b),
            weighted_concordance=weighted_concordance(states_a, states_b),
            shared_vrs_ids=ids_a & ids_b,
            total_a=len(ids_a),
            total_b=len(ids_b),
        )

    def match_against_all(self, context, sample_id, *, top_n=None, candidate_vrs_ids=None):
        if not context.sample_exists(sample_id):
            raise KeyError(sample_id)

        results = [
            self.match_pair(context, sample_id, other, candidate_vrs_ids=candidate_vrs_ids)
            for other in context.list_samples()
            if other != sample_id
        ]
        results.sort(key=lambda result: result.jaccard, reverse=True)
        return results[:top_n] if top_n is not None else results


def create_plugin():
    return MyPlugin()
```

```bash
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --plugin-file my_plugin.py
```

## Development and validation

A good workflow:

1. Start from the example in `examples/plugins/` and adjust only the scoring
   logic first.
2. Run it on a small, known cohort where you understand the expected ranking,
   and compare against the built-in `identity` plugin.
3. Add tests for the new behavior before sharing.

Validate on data where the expected answer is already known. Useful checks:

- known duplicate or resequenced pairs rank near the top;
- clearly unrelated samples rank lower;
- reruns give identical output (keep plugin scoring deterministic);
- `candidate_vrs_ids`-restricted runs behave as expected;
- output still populates a sensible `MatchResult`.

Keep a small regression dataset so you can verify behavior after every change.
Prefer simple, explainable scoring rules over opaque heuristics, especially for
QC or sample-identity decisions, and treat `context` as read-only — never
assume anything about internal DB tables.

## Example plugin ideas

The following are design ideas, not built-in features. They illustrate the
kinds of ranking logic a plugin can implement.

### Rare-variant-weighted identity confirmation

Useful when confirming two samples came from the same donor but you want rare
variants to count for more than common ones — e.g. matching a resequenced
sample to an earlier run, catching a swapped tumor-normal pair, or confirming
identity across sequencing centers.

The default matcher treats all indexed alleles equally, which can understate
rare, high-information variants. Such a plugin could assign each VRS ID a weight
based on cohort frequency, compute a weighted overlap instead of plain Jaccard,
and optionally penalize discordance at rare loci more heavily. Conceptually:

- common allele: weight = 1
- uncommon allele: weight = 3
- rare allele: weight = 10

So samples sharing several rare alleles rank above samples matching mostly on
common ones.

**Data needed:** an external notion of variant frequency — internal cohort
allele frequencies, release-specific site frequencies, panel rarity tiers, or a
hand-maintained allowlist. Weights can be hard-coded for a prototype, loaded
from a sidecar file, or packaged with the plugin.

**Cautions:** rare calls are more sensitive to artifacts; rarity from a small
cohort is unstable; high-confidence loci beat using all rare calls blindly; and
be explicit that `jaccard` now holds a weighted identity score.

### Candidate-gene-only matching

Useful when matching should be driven only by variants in genes, regions, or
panels relevant to a disease or study question — e.g. a cardiomyopathy panel,
inherited-cancer-predisposition genes, or a curated research panel. The question
shifts from genome-wide identity to "are these samples similar *within the genes
I care about*?"

This is often the simplest custom plugin, because the core logic can be
identical to the built-in matcher after filtering. Such a plugin could restrict
scoring to VRS IDs in a supplied gene list (reusing `candidate_vrs_ids` when the
set is already known), load a gene-to-VRS mapping from a sidecar file, and
optionally add rules like weighting loss-of-function candidates more.

**Data needed:** a curated list of candidate VRS IDs, a gene-to-VRS lookup
table, a panel definition, or phenotype-driven regions converted to VRS IDs.

**Cautions:** panel similarity is not genome-wide identity; small candidate sets
are sensitive to missing calls; different panel definitions yield different
rankings; make sure QC users know the matching is panel-restricted.

### Combining the two

A plugin could restrict matching to a candidate gene set *and* weight variants
within it by rarity or predicted impact. This suits phenotype-driven research
but should be validated carefully, since the score becomes more specialized and
less comparable to ordinary genome-wide identity matching.

## Packaging a reusable plugin

To let others install your plugin with `pip`, expose it via Python entry points
in `pyproject.toml`:

```toml
[project.entry-points."vrs_matcher.plugins"]
my_plugin = "my_package.my_plugin:create_plugin"
```

The referenced object may be a plugin instance, a plugin class with a
zero-argument constructor, or a zero-argument factory returning an instance.
After installation, users select it by name:

```bash
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --algorithm my-plugin
```

## Troubleshooting

| Error | Cause |
| --- | --- |
| `Plugin script must define create_plugin()` | Your script file is missing a `create_plugin()` function. |
| `Matcher plugin not found: ...` | The `--algorithm` name doesn't exist, the plugin package isn't installed, or `--algorithm` doesn't match the `name` in a script plugin. |
| `Plugin '...' API version '...' is not supported` | Your plugin's `api_version` doesn't match the current API. Set `api_version = PLUGIN_API_VERSION` (from `vrs_matcher.plugins`). |
