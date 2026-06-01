# Matcher plugins

`vrs-matcher` can run matching algorithms supplied by users, not just the
built-in algorithm. This is useful when you want to:

- try a different similarity score,
- prioritize rare or clinically important variants,
- add lab-specific ranking rules, or
- prototype a method before turning it into a reusable package.

There are three supported plugin sources:

1. **Built-in plugins** shipped with `vrs-matcher`
2. **Entry-point plugins** from installed Python packages
3. **Local script plugins** loaded directly from a `.py` file

For most bioinformatics users, the easiest path is:

1. start with a **local script plugin**,
2. test it on your cohort,
3. then package it as an **entry-point plugin** if you want to share it.

## What a plugin changes

A plugin changes **how samples are scored and ranked after data are already
loaded into the SQLite index**.

Plugins do **not** replace:

- VCF parsing,
- VRS allele extraction,
- SQLite storage, or
- the `MatchResult` output structure.

In other words, plugins operate on indexed sample-level allele/genotype data,
not on raw VCF records.

## When should a bioinformatician write a plugin?

Write a plugin when the built-in identity matcher is close, but not quite the
ranking logic you need.

Common examples:

- **Sample identity confirmation** where rare variants should contribute more
  than common variants
- **Cohort deduplication** where you want to down-weight noisy loci
- **Disease-focused matching** where only a candidate gene list or panel should
  influence ranking
- **QC workflows** where discordant zygosity at high-confidence loci should be
  penalized more strongly than the default matcher does

If your goal is only to restrict the comparison to a subset of variants, first
ask whether `candidate_vrs_ids` is enough. If you need different scoring or
ranking behavior, a plugin is the right tool.

## Quick start

### 1. See what plugins are available

```bash
uv run vrs-matcher plugins list
```

This lists built-in plugins and any installed entry-point plugins.

### 2. Run the built-in algorithm explicitly

```bash
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --algorithm identity
```

### 3. Run a local script plugin

```bash
uv run vrs-matcher match-sample SAMPLE_A \
  --db matches.db \
  --plugin-file examples/plugins/jaccard_floor_plugin.py
```

You can also use plugins with pairwise matching or shared-variant output:

```bash
uv run vrs-matcher match-samples SAMPLE_A SAMPLE_B \
  --db matches.db \
  --plugin-file examples/plugins/jaccard_floor_plugin.py

uv run vrs-matcher shared-variants SAMPLE_A SAMPLE_B \
  --db matches.db \
  --plugin-file examples/plugins/jaccard_floor_plugin.py
```

## How plugin selection works

### Built-in or installed plugin

Use `--algorithm`:

```bash
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --algorithm identity
```

### Local script plugin

Use `--plugin-file`:

```bash
uv run vrs-matcher match-sample SAMPLE_A \
  --db matches.db \
  --plugin-file my_plugin.py
```

When using `--plugin-file`, you usually do **not** need `--algorithm`.

If you do provide `--algorithm`, it must match the plugin's `name`.

## Plugin contract

A plugin object must define all of the following:

- `name`: unique plugin name, for example `"identity"` or `"rare-priority"`
- `api_version`: currently `"1"`
- `match_pair(context, sample_a, sample_b, *, candidate_vrs_ids=None) -> MatchResult`
- `match_against_all(context, sample_id, *, top_n=None, candidate_vrs_ids=None) -> list[MatchResult]`

For **script plugins**, the Python file must also define:

- `create_plugin()`: returns the plugin object

### What the plugin receives

Plugins are given a read-only `context` object instead of direct SQL access.
This keeps plugin code simple and decoupled from the database schema.

This means your plugin can assume that:

- samples have already been loaded,
- VRS IDs have already been normalized,
- genotype states are already available per sample, and
- the plugin only needs to focus on scoring logic.

Available helpers:

- `context.sample_exists(sample_id)`
- `context.get_vrs_ids(sample_id)`
- `context.get_genotype_states(sample_id)`
- `context.list_samples()`

### What the plugin returns

Plugins return `MatchResult` objects from `vrs_matcher.matcher`.

That means your plugin should populate:

- `sample_a`
- `sample_b`
- `jaccard`
- `weighted_concordance`
- `shared_vrs_ids`
- `total_a`
- `total_b`

Even if your algorithm is not truly Jaccard-based, the current CLI and notebook
interfaces still expect the score to be reported in the `jaccard` field. If you
use that field for a custom primary score, document that clearly for users.

## Recommended development workflow

1. Start from the example plugin in `examples/plugins/`
2. Rename the plugin and adjust only the scoring logic first
3. Run it on a small, known cohort where you already understand the expected ranking
4. Compare results against the built-in `identity` plugin
5. Add tests for the new behavior before sharing it with others

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

Run it with:

```bash
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --plugin-file my_plugin.py
```

## How to validate a new plugin

Before using a plugin in production or in a release-QC workflow, validate it on
data where the expected answer is already known.

Recommended checks:

- a known duplicate or resequenced sample pair should rank near the top
- clearly unrelated samples should rank lower
- rerunning the same command should give the same output
- candidate-gene or panel-restricted runs should behave as expected when
  `candidate_vrs_ids` is supplied
- plugin output should still populate a sensible `MatchResult`

Useful commands while iterating:

```bash
uv run vrs-matcher plugins list

uv run vrs-matcher match-samples SAMPLE_A SAMPLE_B \
  --db matches.db \
  --plugin-file my_plugin.py

uv run vrs-matcher match-sample SAMPLE_A \
  --db matches.db \
  --top 5 \
  --plugin-file my_plugin.py
```

## Worked example

The repository includes a complete example plugin:

- `examples/plugins/jaccard_floor_plugin.py`

That plugin behaves like the built-in identity matcher, but it enforces a
minimum Jaccard score floor. It is meant as a simple teaching example, not
necessarily as a recommended production ranking method.

Run it with:

```bash
uv run vrs-matcher match-sample SAMPLE_A \
  --db matches.db \
  --plugin-file examples/plugins/jaccard_floor_plugin.py
```

## Example use-case ideas for custom plugins

The examples below are intended as design ideas for custom plugins. They are
not built into `vrs-matcher`, but they illustrate the kinds of ranking logic a
plugin can implement.

### Rare-variant-weighted identity confirmation

#### When this is useful

This idea is useful when you are trying to confirm that two samples came from
the same donor, but you want **rare variants to contribute more evidence than
common variants**.

Examples:

- confirming that a resequenced sample matches an earlier run,
- checking whether a tumor-normal pair was accidentally swapped,
- confirming sample identity across sequencing centers or pipelines,
- strengthening identity comparisons in cohorts where many common variants are
  shared across individuals.

#### Why a rare-variant-weighted plugin can help

The default identity matcher treats all indexed VRS alleles equally in its set
overlap calculation. In some workflows, that may understate the value of rare,
high-information variants.

For example:

- sharing a very common variant may provide only weak evidence of identity,
- sharing a rare singleton or near-singleton variant may provide much stronger evidence,
- discordance at a rare, high-confidence site may be more concerning than
  discordance at a common site.

#### What the plugin might do

A rare-variant-weighted plugin could:

- assign each VRS ID a weight based on cohort frequency,
- compute a weighted overlap score instead of plain Jaccard,
- optionally increase the penalty for genotype discordance at rare loci,
- still return a standard `MatchResult`, using the `jaccard` field to store the
  custom primary score.

Conceptually, you might score alleles like this:

- common allele: weight = 1
- uncommon allele: weight = 3
- rare allele: weight = 10

Then two samples sharing several rare alleles would rank above two samples that
match mostly on common alleles.

#### Data needed

This kind of plugin usually needs an external or precomputed notion of variant
frequency, for example:

- internal cohort allele frequencies,
- release-specific site frequencies,
- panel-specific rarity tiers, or
- a simple hand-maintained allowlist of highly informative loci.

Those weights could be hard-coded for a prototype, loaded from a sidecar file,
or packaged with the plugin.

#### Things to be careful about

- Rare variant calls may be more sensitive to calling artifacts.
- If rarity is computed from a small cohort, weights can be unstable.
- Highly filtered, high-confidence loci are usually better than using all rare calls blindly.
- Be explicit that the returned `jaccard` field is now a weighted identity score,
  not a literal Jaccard coefficient.

#### Good validation questions

- Do known duplicate or resequenced samples move higher in the ranking?
- Do unrelated samples remain well separated?
- Are rankings robust if a few noisy rare calls are removed?
- Does the plugin behave sensibly across different cohorts or sequencing batches?

### Candidate-gene-only sample matching

#### When this is useful

This idea is useful when the matching task should be driven only by variants in
genes, regions, or panels relevant to a particular disease or study question.

Examples:

- restricting comparison to a cardiomyopathy gene panel,
- comparing samples only within inherited cancer predisposition genes,
- prioritizing matching based on phenotype-relevant genes,
- checking whether two samples share the same profile within a curated research panel.

#### Why a candidate-gene-only plugin can help

Sometimes genome-wide identity is not the only question. You may want to ask:

- are these samples similar **within the genes I care about**?
- do they share the same signal in a phenotype-driven subset of loci?
- does a ranking change if I ignore the rest of the genome?

This can be especially useful in exploratory or disease-focused research
settings, where phenotype relevance matters more than whole-genome similarity.

#### What the plugin might do

A candidate-gene plugin could:

- restrict scoring to VRS IDs mapped to a supplied gene list,
- use the existing `candidate_vrs_ids` input if the candidate set is already known,
- load a gene-to-VRS mapping from a sidecar file,
- compute the same score as the default matcher, but only over the candidate subset,
- optionally rank by additional rules, such as giving loss-of-function candidate variants more weight.

In practice, this is often the simplest custom plugin to build because the core
logic may be identical to the built-in matcher after filtering to a subset.

#### Data needed

You typically need one of the following:

- a curated list of candidate VRS IDs,
- a gene-to-VRS lookup table,
- a panel definition generated upstream,
- or phenotype-driven region definitions converted to VRS IDs before matching.

#### Things to be careful about

- Similarity within a small panel is not the same as genome-wide identity.
- A small candidate set can make scores more sensitive to missing calls.
- Different panel definitions can produce very different rankings.
- If you use the plugin for QC, make sure users understand it is panel-restricted matching.

#### Good validation questions

- Do samples with known phenotype-relevant overlap rise in rank as expected?
- Does the result remain interpretable when only a small number of loci are used?
- Are missing or low-quality candidate loci affecting the ranking too strongly?
- Is the plugin still useful when the candidate panel changes between projects?

### Combining the two ideas

These two ideas can also be combined.

For example, a plugin could:

- restrict matching to a candidate gene set, and then
- weight variants within that set by rarity or predicted impact.

That kind of plugin may be useful for phenotype-driven research workflows, but
it should be validated carefully because the score becomes more specialized and
less comparable to ordinary genome-wide identity matching.

## Packaging a reusable plugin

If you want others to install your plugin with `pip`, expose it via Python
entry points in `pyproject.toml`:

```toml
[project.entry-points."vrs_matcher.plugins"]
my_plugin = "my_package.my_plugin:create_plugin"
```

The referenced object may be:

- a plugin instance,
- a plugin class with a zero-argument constructor, or
- a zero-argument factory returning a plugin instance.

After installation, users can call it by name:

```bash
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --algorithm my-plugin
```

## Practical recommendations

- Start by editing a local script plugin and validating it on a small cohort.
- Keep plugin code deterministic: the same inputs should produce the same scores.
- Treat `context` as read-only; do not assume anything about internal DB tables.
- Reuse `candidate_vrs_ids` if your workflow includes candidate-gene or panel-restricted matching.
- If your algorithm changes ranking semantics, document that clearly for downstream users.
- Prefer simple, explainable scoring rules over opaque heuristics when the plugin will be used for QC or sample identity decisions.
- Keep a small regression dataset so you can verify behavior after every plugin change.

## Troubleshooting

### `Plugin script must define create_plugin()`

Your script file is missing a `create_plugin()` function.

### `Matcher plugin not found: ...`

Either:

- the `--algorithm` name does not exist,
- the plugin package is not installed, or
- `--algorithm` does not match the `name` inside a script plugin.

### `Plugin '...' API version '...' is not supported`

Your plugin's `api_version` does not match the current plugin API expected by
`vrs-matcher`.

Set:

```python
api_version = PLUGIN_API_VERSION
```

