# ADR: Benchmark vrs-matcher against VCFtools for sample comparison

## Status

Proposed · 2026-09-30. This document specifies future work; no benchmark has been run and no superiority claim is made.

Repository baseline reviewed: `9382c731bad1335fcdc47c31e71fd44353bc30b6`.

Implementation update: a separate [KING relatedness track](user-story-king-robust-plugin.md)
now has an implementation and [validation report](king-validation.md). Its
`scripts/benchmark_king.py` runner compares against `--relatedness2`; it does not
implement or replace this ADR's proposed genotype-discordance experiment.

## Context

An analyst asked how this repository compares to VCFtools. The useful overlapping task is sample comparison, not the full breadth of VCF processing. The repo's default `identity` plugin compares carried VRS allele sets and ranks indexed samples. VCFtools individual discordance provides a genotype-based baseline for supplied sample correspondences. Their native scores have different denominators; equating weighted concordance to one minus VCFtools discordance would be incorrect.

See the [bioinformatics comparison](comparison-vcftools.md), [current matching implementation](../src/vrs_matcher/matcher.py), and [VCFtools individual-discordance implementation](https://github.com/vcftools/vcftools/blob/master/src/cpp/variant_file_diff.cpp).

## Decision

Use **C++ `vcftools --diff-indv-discordance`** as the primary baseline. Pin the actual executable version and source revision; the published manual describes v0.1.16. Explicitly exclude the separate Perl `vcf-compare` utility and bcftools from this first comparison. [VCFtools manual](https://vcftools.github.io/man_latest.html)

Evaluate three independent questions:

1. **Correctness:** Does each tool implement its own stated score on controlled inputs?
2. **Utility:** How reliably does each score recover labeled sample identity under defined perturbations?
3. **Cost:** What are preparation, first-result, and repeated-query time, memory, and storage costs?

Start with one opt-in benchmark script and small fixtures, using the existing Python, SQLite, pysam, and VRS integration dependencies. Keep VCFtools an external benchmark dependency. Do not change production scoring, add a framework, or implement a new plugin to force score equality.

## Input contract and fair comparison

The primary comparison uses one reference assembly, autosomal diploid biallelic SNVs, identical sorted loci and REF/ALT strings, and fully called genotypes. Materialize a common site panel and sample selection before either tool runs. Keep `0/0` calls for VCFtools: removing them to imitate the index would change the baseline's meaning.

Use a common preparation pass for record filtering and genotype QC. Retain PASS/unfiltered records; mask calls with present GQ <20 or DP <0. Missing QC fields are allowed, matching the loader. For the primary panel, retain only loci called in every selected observation after masking, so no pair benefits from a different callable denominator. Save excluded-record and excluded-call counts. Do not assume similarly named tool flags have identical missing-value semantics. This policy follows [loader.py](../src/vrs_matcher/loader.py).

For VRS inputs, verify exactly one valid ALT ID per retained biallelic record, with no missing/placeholder IDs. For later multiallelic tests, verify ALT-to-ID order and cardinality. Fail preparation on annotation errors instead of silently comparing a reduced subset. Record annotation failures separately in robustness experiments.

Assign distinct observation IDs such as `A__HG00096` and `B__HG00096` before indexing; store the true donor ID separately. The database key is `(sample_id, vrs_id)`, not `(source_dataset, sample_id, vrs_id)`. Create a fresh database for every build measurement. Verify actual SQL row counts, rather than interpreting the loader's emitted-row count as the final unique-row count. [Database implementation](../src/vrs_matcher/db.py)

## Data and experimental tracks

### 1. Small, hand-checkable correctness fixtures

Author a compact VCF fixture with a recorded truth table, independent of either scoring implementation. Cover identical nonempty samples, disjoint ALT sets, `0/0` agreement, `0/1` versus `1/1`, reordered phased/unphased diploid calls, fully missing calls, `./1`, multiallelic `1/2`, failed FILTER, low/missing GQ/DP, unannotated records, and empty samples. Run exceptional cases separately from the primary panel.

Assert native scores and counts against the hand-computed expectations. Include the three-locus example in the comparison document: Jaccard 1, weighted concordance 0.75, and VCFtools discordance 1/3. For VCFtools partial-call and multiallelic behavior, inspect the pinned version and record its actual comparison denominator; do not extrapolate from the complete diploid case. [VCFtools source](https://github.com/vcftools/vcftools/blob/master/src/cpp/variant_file_diff.cpp)

Treat a pair with zero common called loci as unscorable in the baseline. Record the matcher's empty-set Jaccard of 1, but classify a zero-evidence identity decision as **insufficient evidence** in the benchmark wrapper. Preserve raw output so this wrapper policy cannot hide implementation behavior.

### 2. Public-data identity experiment

Reuse the source URL and 20 AFR/EUR sample IDs in [the existing fixture](../tests/integration/conftest.py). Its fetch interval is chr22 `[16,000,000, 17,000,000)` in pysam's zero-based coordinates. Freeze the downloaded slice locally and hash it; reuse fixture logic without importing or invoking pytest fixtures from the benchmark runner.

Create a reference cohort and separately named query copies. Exact copies are positive controls. Add deterministic, labeled perturbations using a recorded seed: 0%, 1%, 5%, and 10% missing calls or genotype substitutions, in separate experiments. At a selected biallelic site, substitute uniformly among the other two diploid genotypes. Include unrelated-donor queries with no true match in the reference cohort. These are synthetic replicates of public genotypes, not independent sequencing replicates or a contamination model.

Select calibration and evaluation donors before making copies; keep every copy of a donor in the same partition. Include same-population nonmatches and, where independently documented, relatives as difficult negatives. Do not assume every distinct 1000 Genomes sample is unrelated. Keep the existing population-structure assertion as a separate sanity check, not the identity endpoint. Broad accuracy claims require a subsequent independently generated replicate-callset dataset with verified donor labels.

### 3. Representation and missingness stress tests

After the primary track passes, add small indels with verified equivalent encodings, ALT reordering, split/unsplit multiallelic records, partial calls, and unequal assay footprints. Verify equivalent edited sequences independently before labeling a representation-only transformation as preserving an allele. Generate VRS annotations for each representation independently; do not copy IDs to manufacture agreement.

For representation experiments, compare both original encodings and a harmonized representation supplied to both tools; pin the normalizer and reference and account for that preparation time. Unsupported transformations and annotation failures are explicit outcomes. Cross-assembly, symbolic SV, and haplotype-equivalence comparisons are deferred.

For missingness experiments, report two conditions: an oracle common-callable panel derived from both observations, and the unharmonized inputs. Label the former as requiring extra upstream information. The SQLite index cannot reconstruct callable reference loci or distinguish missing calls from absent alleles. Do not interpret reduced Jaccard automatically as biological disagreement.

## Execution protocol

Prepare two cohort files with unique observation names and a donor-label manifest. For a known pair, create a two-column, tab-delimited `pair.tsv` containing the name in file A followed by the name in file B. A representative baseline invocation is:

```bash
vcftools --gzvcf prepared-a.vcf.gz --gzdiff prepared-b.vcf.gz \
  --diff-indv-map pair.tsv --diff-indv-discordance --out results/pair
```

The runner must select the mapped row in `.diff.indv`, verify its `N_COMMON_CALLED` and `N_DISCORD`, and preserve logs. Directory and input preparation are prerequisites for this command. Map orientation and output naming follow the [VCFtools manual](https://vcftools.github.io/man_latest.html).

For cohort search, compare every query against the same candidate IDs using repeated one-to-one mappings; exclude unmapped rows. Batch independent pairs in a mapping where supported, and measure this batched known-pair workload separately. Do not treat `--diff-discordance-matrix` as an all-samples similarity matrix. The search wrapper and its file-processing costs belong in baseline search timing.

Load independently annotated versions into a fresh index:

```bash
uv run vrs-matcher load-samples prepared-a.vrs.vcf.gz --db results/run.db
uv run vrs-matcher load-samples prepared-b.vrs.vcf.gz --db results/run.db
uv run vrs-matcher match-samples A__HG00096 B__HG00096 --db results/run.db
```

For scoring accuracy, use the Python API to save full-precision values; CLI scores are rounded. For search timing, load the candidate cohort plus one query and invoke `match_against_all(..., top_n=None)`, so other queries do not accidentally become candidates. Record query annotation/loading cost separately. Measure one query versus N candidates and all unique pairs as distinct workloads. `top_n` truncates after all peers are scored; it does not bound the current search work. [Matcher](../src/vrs_matcher/matcher.py), [CLI](../src/vrs_matcher/cli.py)

## Metrics and interpretation

Save pair-level counts: carried IDs in each sample, intersection, union, shared zygosity agreements/disagreements, common called loci, discordant loci, excluded calls, and excluded records. Preserve native Jaccard, weighted concordance, and `D = N_DISCORD / N_COMMON_CALLED`; use `1-D` only as a clearly labeled baseline concordance score.

On the common biallelic panel, calculate a simple independent genotype contingency table from prepared VCFs. It should reproduce each native score under its own definition. This checks denominator accounting without making either program the truth oracle.

For identity retrieval, rank by Jaccard for the built-in matcher and by `1-D` for VCFtools. Report top-1 accuracy, recall@5, tie frequency, and the gap to the best incorrect candidate. Do not count arbitrary ordering within a tie as successful unique identification. Report weighted concordance as a secondary score; any combined decision rule must be fixed using calibration data only.

For same-donor decisions, calibrate separate thresholds and minimum-evidence counts for each score, freeze them, then report sensitivity, false-positive rate, precision, and abstention rate on evaluation donors. Include pair and donor counts and donor-resampled uncertainty intervals; pairs sharing a donor are not independent. Stratify by perturbation, overlap, and population. With 20 donors, results are a pilot, not evidence for a production false-match rate.

## Performance and reproducibility

Run on one documented machine, sequentially, with the same local inputs and resource limits. Record CPU, RAM, OS, filesystem/storage, process/thread settings, repository SHA, lockfile hash, Python/SQLite/cyvcf2/pysam/VRS versions, VCFtools build, SeqRepo snapshot/reference digests, input hashes, preparation commands, sample lists, seed, and parameter manifest. Downloads and environment installation are setup costs, recorded separately from offline compute.

Time these stages separately:

| Stage | vrs-matcher | VCFtools baseline |
|---|---|---|
| Shared preparation | Subsetting, QC, harmonization | Same prepared observations |
| Tool-specific preparation | VRS annotation and SQLite build | Mapping and any materialized comparison files |
| First result | Preparation plus one comparison/search | Preparation plus one comparison/search |
| Repeated work | Reuse index; include new-query ingestion when applicable | Reuse prepared files; include repeated scans and wrapper overhead |

Use a fresh process for each timed run. Record one first-run observation and five subsequent repetitions, alternating tool order; report all measurements plus median and range. Fresh process does not mean cold filesystem cache: label cache conditions honestly, and do not claim cold-cache results without controlling the cache. Use platform `time` resource reporting or equivalent subprocess resource accounting for wall time, CPU time, and peak RSS, explicitly converting platform-specific RSS units. Record database, annotated VCF, mapping, and temporary-file sizes separately from shared reference storage.

Start with the 20-sample slice. Then, as capacity permits, use fixed subsets of 100 and 500 donors and nested 10k/100k/1M-SNV panels from the same pinned source, recording actual available counts. Vary samples and variants separately. Record timeouts, memory failures, annotation failure rates, and incomplete runs; never omit unsuccessful configurations from the report.

For Q searches, compare measured totals `T_shared + T_annotation + T_index + sum(T_query_load + T_search)` against `T_shared + T_baseline_prepare + sum(T_baseline_search)`. Report any observed amortization crossover and the tested range; do not infer one from a faster warm query alone.

## Deliverables and completion criteria

The future implementation should produce a manifest, deterministic fixture inputs and expected counts, pair-level score TSV, timing/resource TSV, raw command logs, and a short Markdown report. Store large public inputs and SQLite artifacts outside version control. Add one runnable correctness check using the existing test setup; keep performance runs opt-in and out of normal CI timing assertions.

The benchmark is complete when both tools pass the hand-checked cases under their documented semantics, pair counts and candidate lists reconcile, all planned pilot cases have results or explicit failure reasons, and another engineer can reproduce the report from pinned inputs. Numeric equality between unlike scores is not an acceptance condition. Any accuracy or speed claim must identify its workload, evidence denominator, and included preparation costs.

## Alternatives and consequences

- **Compare native score numbers directly:** rejected because shared-allele agreement and called-site genotype agreement answer different questions.
- **Use only the existing population test:** rejected as an identity benchmark; population similarity can make nonmatches appear close.
- **Time only an already-built index:** retain as a repeated-query measurement, but it cannot establish end-to-end superiority.
- **Build a general benchmarking platform or tune a new matcher:** deferred until this small experiment identifies a concrete need.

This decision yields an interpretable comparison with modest implementation work. Its limits remain explicit: synthetic errors cannot reproduce all assay effects, pilot cohort size limits accuracy estimates, and the VCFtools search wrapper is part of the measured workflow rather than a native cohort-search feature.
