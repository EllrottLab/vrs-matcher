# Estimate kinship with the KING-robust plugin

`king-robust` estimates relationships using explicit, jointly called genotypes
on a declared human autosomal SNP panel. It reports the Manichaikul
between-family estimate (equation 11) for ranking, the within-family estimate
(equation 9) as a diagnostic, and the supporting counts. This is an estimator,
not an automatic pedigree or clinical identity classifier.

## Prepare the inputs

Use one reference assembly and consistently annotated, diploid biallelic SNVs.
VRS generation and sequence normalization happen upstream. The index checks
identifier syntax and agreement with the declared panel, but does not recompute
VRS IDs or independently verify your reference sequence. An identifier must be
`ga4gh:VA.` followed by a 32-character base64url digest; reference sequence IDs
use `ga4gh:SQ.` with the same digest length.

Create a tab-delimited panel with these exact column names (column order is
flexible). Positions are **one-based VCF positions**:

```text
reference  annotation_version  chrom  sequence_id  pos  ref  alt  vrs_id
```

Every row must contain all eight values. `reference` and `annotation_version`
are consistent, nonempty provenance identifiers chosen for this panel; record
the actual assembly/snapshot and annotation software version. `chrom` is a human
autosome (`1` through `22`, optionally prefixed with `chr`). `sequence_id` is its
reference sequence digest. `ref` and `alt` are distinct uppercase A/C/G/T bases.
Exactly one marker is allowed at a reference sequence/position, and every ALT
VRS ID must be unique. One genotype index uses one panel and one QC policy.
The panel's SHA-256 file hash is its identifier; reuse the same file for appends.

Declare matching provenance and sequence identities in the annotated VCF header:

```text
##reference=YOUR_REFERENCE_IDENTIFIER
##vrs_annotation=YOUR_ANNOTATION_VERSION
##contig=<ID=chr1,length=YOUR_LENGTH,refget=YOUR_GA4GH_SQ_IDENTIFIER>
##INFO=<ID=VRS_Allele_IDs,Number=A,Type=String,Description="VRS ALT allele IDs">
```

`YOUR_*` values above are placeholders. `1` and `chr1` are accepted as aliases
only when the used VCF contig's `refget` identity agrees with the panel. Other
contig aliases must be harmonized upstream. Both header provenance values must
match the panel exactly. Each eligible panel record must contain exactly its
expected ALT ID in `INFO/VRS_Allele_IDs`.

The [small panel](../tests/data/king/panel.tsv) and
[VCF](../tests/data/king/oracle.vcf) are runnable **synthetic fixtures** with
synthetic identifiers, not biological reference data.

## Index and query

Use a new database and unique observation names. Independently measured or
reprocessed observations of the same donor must have different names:

```bash
uv run vrs-matcher load-samples cohort.vrs.vcf.gz --db cohort.db \
  --index-genotypes --panel king-panel.tsv --source-dataset release-1
uv run vrs-matcher match-samples OBSERVATION_A OBSERVATION_B \
  --db cohort.db --algorithm king-robust --json
uv run vrs-matcher match-sample OBSERVATION_A --db cohort.db \
  --algorithm king-robust --top 10 --json
```

For a local demonstration, replace the first input with
`tests/data/king/oracle.vcf`, the panel with `tests/data/king/panel.tsv`, and
observation names with `A` and `B`. Expected kinship is 0.0 and the within-family
diagnostic is 0.1. The synthetic example is intentionally too small for
biological inference.

Ingestion returns a JSON summary and writes both indexes in one transaction.
Failed loads roll back all new observations and calls. Successful observations
are immutable: appending an existing name is rejected. To replace a sample or
change panel/QC policy, build a new database from source VCFs.

The genotype table stores passing dosages **0, 1, and 2**. An absent row means
unknown/excluded, never reference. Fully missing, partially missing, invalid,
and non-diploid genotypes are excluded from KING. Indels, multiallelic records,
non-panel records, and non-panel contigs are excluded and counted. Duplicate
panel records, inconsistent marker orientation, and missing/invalid annotations
at eligible panel records fail the load.

Default QC retains PASS/unfiltered records, GQ ≥20 and DP ≥0 when those fields
exist. Missing QC is allowed and counted. `--gq` and `--dp` apply equally to
reference and alternate calls. The new path's per-record/call summaries are in
`genotype_observation.summary`; provenance and input hashes are in
`genotype_observation.provenance`. Record exclusions are counted once per VCF,
call exclusions per sample; missing QC counters overlap call outcomes.

The carried-allele index still follows the existing identity projection, including
its partial-call handling. KING-specific exclusions do not redefine identity.
Reference genotypes never become carried VRS IDs. Existing allele-only databases
still support `identity`, but **cannot be upgraded by inference**: re-ingest the
original inputs into a new genotype-enabled index. A KING query involving an
allele-only observation fails with a re-ingestion message.

## Read the results

Pair results expose full precision with `--json`:

| Field | Interpretation |
|---|---|
| `panel_id`, `estimator`, `estimator_version` | Reproducibility identifiers |
| `n_common` | M: jointly called panel markers |
| `het_a`, `het_b` | H_i, H_j on the jointly called markers |
| `het_both` | HH: both samples heterozygous |
| `opposite_hom` | O: opposite homozygotes, or IBS0 count |
| `kinship` | `0.5 - (H_i + H_j - 2*HH + 4*O)/(4*min(H_i,H_j))` |
| `kinship_within_family` | `(HH - 2*O)/(H_i + H_j)` |
| `ibs0_fraction` | O/M |
| `status`, `reason`, `within_family_reason` | Explicit explanation of unscorable estimates |

Negative estimates remain negative. Zero denominators yield JSON `null` (text
`NA`), not zero kinship. The two estimates can have different denominator
availability. Empty observations never yield a perfect relationship score.

Search returns `KinshipMatches` with `matches` and `unscorable` lists; JSON uses
the same keys. Only finite between-family estimates enter top-N, ordered by
score descending then sample ID. Unscorable peers are reported separately and
are not discarded by top-N. Pair queries allow self-comparison; cohort search
excludes the query itself. `shared-variants` rejects kinship results.

Python APIs accept `candidate_vrs_ids` to restrict comparisons to a panel
subset, including its explicit reference calls. Unknown IDs are errors and an
empty restriction is unscorable. The result types live in
`vrs_matcher.models`. Existing `identity` APIs still return `MatchResult` and
lists of `MatchResult`.

## Validate and benchmark

Run deterministic acceptance tests with the normal suite:

```bash
uv run pytest tests/test_king.py
```

For the opt-in external comparison, supply a VCFtools **0.1.16 or newer** executable:

```bash
VCFTOOLS=/path/to/vcftools uv run pytest \
  tests/integration/test_king_vcftools.py --run-integration
```

This tests exact underlying counts and within-family numerical agreement on
complete calls, plus the documented missingness divergence. It does not demand
that VCFtools equal the plugin's between-family estimate.

The benchmark runner requires Linux or macOS (`wait4` resource accounting),
VCFtools, and the project environment. It uses fresh processes, six repetitions,
alternating load/baseline order, and an uncontrolled OS cache:

```bash
uv run python scripts/benchmark_king.py \
  --vcf tests/data/king/oracle.vcf --panel tests/data/king/panel.tsv \
  --vcftools /path/to/vcftools --output /tmp/king-smoke --smoke
```

The output directory must not exist. Inputs must already be restricted to panel
markers and QC-masked; low GQ/DP calls and failed records are rejected rather
than giving VCFtools a different input. Complete diploid calls or `./.` are
accepted. Missing-call results are retained as estimator differences.

Outputs include manifests with input/source/executable hashes, pair estimates,
VCFtools comparisons, raw logs, timing/RSS TSV, database sizes, and a Markdown
report. The report separates index construction (`load`),
one-versus-all matching (`query`), all distinct unordered pairs (`all`), and
VCFtools `--relatedness2`. Separate relative percentages compare indexed
all-pairs matching and index-plus-all-pairs against the direct VCFtools run.
These are different workloads: indexed matching excludes VCF reading for
VRS-Matcher, while the direct VCFtools run includes it; index-plus-all-pairs
includes VRS-Matcher indexing. The report explains these caveats. It shows
repeat 0 separately, then the median and range of later runs; peak RSS is the
maximum across runs for each stage. See the
[performance ADR](adr-king-performance.md) for planned architecture
optimizations and their validation criteria.

The report includes separate relative-performance percentages for indexed
all-pairs matching and index-plus-all-pairs: `100 × (VCFtools median time /
VRS-Matcher median time)`, with VCFtools as the 100% reference. Below 100%
means VRS-Matcher took longer for that reported workload; above 100% means it
was faster. The indexed comparison excludes VCF reading for VRS-Matcher, while
the VCFtools run includes reading the prepared VCF; the pipeline comparison
includes VRS-Matcher index construction. Neither is a pure algorithm-only
comparison. VCFtools computes ordered pairs and diagonals, whereas the Python
all-pairs worker computes each distinct unordered pair once. The query worker
measures one-versus-all separately. All measurements include process startup
but **exclude upstream preparation, VRS annotation, downloads, and environment
setup**; record those costs separately before making an end-to-end claim. The
runner stores all pair results in memory; large-cohort runs may need streamed
output after profiling establishes that need.

For biological evaluation, omit `--smoke` and supply both:

- `--samples samples.tsv`: columns `sample`, `donor`, `family`, `ancestry`, one
  row per indexed observation. Use a distinct family block for each unrelated
  donor and the same block for known relatives/replicates.
- `--relationships relationships.tsv`: columns `sample_a`, `sample_b`,
  `relationship`, one row per independently verified pair. The evaluation
  must include `duplicate`, `parent_offspring`, `siblings`, and `unrelated`.

The runner produces separate summaries for both plugin estimators and the actual
VCFtools phi, stratified by relationship, ancestry, and callable overlap, with common-locus ranges,
duplicate retrieval with ties reported, and family-block bootstrap intervals
for retrieval accuracy (seed 0; unavailable with fewer than two blocks). No
threshold is trained. Any future calibration must keep families disjoint from
evaluation. Label verification, independent replicates, representative ancestry
coverage, genome-wide marker selection, and annotation provenance remain the
analyst's responsibility. Numerical agreement with VCFtools is not biological
validation. See the [validation report](king-validation.md) for work actually run
and the [user story](user-story-king-robust-plugin.md) for acceptance criteria.
