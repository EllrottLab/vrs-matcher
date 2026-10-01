# ADR: GREGoR phenotype snapshots with batched result enrichment

Status: Proposed · 2026-09-30. This document specifies implementation and
validation; it does not claim the feature or performance targets are achieved.

## Context

The [phenotype user story](gregor-phenotype-use-case.md) requires participant
phenotypes to support case comparison, family review, metadata reconciliation,
cohort selection, and reanalysis. The current [SQLite index](index-description.md)
identifies genomic observations by `sample_id`; it has no participant mapping.
Multiple observations may represent one participant.

Phenotypes change independently of genotype calls. They must retain explicit
positive, negative, and unknown assertions, as well as missing data and source
conflicts. They must not alter identity scores, KING estimates, or allele keys.
The [KING performance ADR](adr-king-performance.md) establishes batched retrieval,
tiled all-pairs scoring, and bounded top-N selection that this work must preserve.

## Decision

Use an optional **SQLite phenotype sidecar**, containing immutable import
snapshots, explicit observation mappings, and participant-level assertions.
Enrich results after genetic scoring and selection. Use standard-library SQLite,
CSV, JSON, and hashing; introduce no service, ORM, ontology engine, or new matcher
plugin. Keep the existing plugin protocol and genetic result dataclasses intact.

A genotype-only command must neither open the sidecar nor execute phenotype SQL.
Phenotype ingestion must not write to the genetic database or require VCF
re-ingestion. Sidecar readers open SQLite in read-only mode; only an explicit
phenotype import command creates or changes its schema.

### Snapshot and linkage model

One snapshot is a complete declared analysis bundle, possibly containing several
source datasets. Queries select exactly one snapshot; there is no implicit
"latest" selection or union of historical snapshots. Each source has a stable,
user-supplied namespace independent of its filename.

Use these tables; all key components are NOT NULL:

| Table | Key and content |
|---|---|
| `phenotype_snapshot` | PK `snapshot_id`; canonical manifest JSON, import timestamp, importer version. Manifest contains GREGoR model version and schema digest, source namespaces, hashes of all input files including mappings/reconciliation, optional ontology release identifiers, and genetic-index binding described below. |
| `phenotype_participant` | Composite PK `(snapshot_id, source_id, participant_id)`; `person_key`, optional source-scoped family ID, and original participant context as JSON. FK to snapshot. |
| `observation_participant` | PK `(snapshot_id, sample_id)`; `source_id`, `participant_id`; composite FK to participant. Many observations can reference one participant. |
| `phenotype_assertion` | PK `(snapshot_id, assertion_id)`; source and participant keys, ontology, term ID, presence, optional source phenotype ID, original fields as JSON, source-file digest, and row number. Composite FK to participant. CHECK presence in `Present`, `Absent`, `Unknown`. |

`person_key` defaults to a collision-safe encoding of the source/participant
pair. An optional explicit reconciliation file assigns the same person key to
cross-source records known to represent one person. Validate that every scoped
participant has exactly one such assignment. Mapping an observation to one
participant anchor then permits retrieval of assertions for all records with
that person key. Preserve their separate sources; reconciliation must not merge
or overwrite assertions. Family identifiers remain source-scoped unless an
explicit family reconciliation is added in a later change.

Create indexes for:

- Participant lookup: `(snapshot_id, person_key)`.
- Assertion retrieval: `(snapshot_id, source_id, participant_id)`.
- Exact cohort selection: `(snapshot_id, ontology, term_id, presence, source_id, participant_id)`.
- Reverse observation lookup: `(snapshot_id, source_id, participant_id, sample_id)`.

SQLite enforces in-sidecar foreign keys on importer connections. There is no
cross-file FK to `samples`; validate mapped observation names against a read-only
view of the genetic registry during import. Do not drop or silently guess
unresolved mappings. Unmapped indexed observations remain valid and report as
unmapped; duplicate/conflicting mapping assignments fail the entire import.

Bind a snapshot to an externally supplied immutable genetic-index artifact
identifier and SHA-256 digest, recording its sample-registry digest as well.
A verified manifest accompanies the genetic artifact. Validate the full file
hash during import or explicit verification, not every query. Query reports
record the binding and actual genetic run provenance; callers must select the
matching verified artifact. A cheap registry comparison can detect name drift
but must not be presented as proof of unchanged genotype content. Replaced or
extended genetic artifacts require explicit rebinding in a new snapshot.
Never rely on a filesystem path as the identity of the genetic data.

### Acquiring GREGoR phenotype and VCF inputs

Acquisition is an explicit preparation step, outside the importer, matcher, and
normal test run. The initial implementation documents an operator-driven export
and download procedure using Terra and the cloud CLI; it does not add a custom
cloud client to the Python package. No data have been downloaded as part of this
ADR. Actual workspace and object identifiers must come from authorized access.

GREGoR participant-level data are controlled-access through AnVIL, under dbGaP
study `phs003047`. Pin a release and its consent group rather than resolving
"latest" at runtime. The release page reviewed on 2026-09-30 lists R05 with
model 1.11, whereas the user story reviewed the development model 1.12. These
versions must not be substituted silently. The release includes a chromosome-
partitioned joint callset, but not every phenotyped participant necessarily has
a sample in it. [GREGoR release and access guidance](https://gregorconsortium.org/data).

The operator obtains the applicable data access approval, links the authorized
AnVIL/Terra account to eRA Commons, and configures a billing project and analysis
workspace. GREGoR workspaces are organized by consent group and use requester
pays. Stage data within the approved analysis environment; downloading to a
workstation is appropriate only when the applicable authorization permits it.
[GREGoR AnVIL setup](https://gregorconsortium.org/research-home/anvil-resources).

#### Phenotype and linkage tables

1. Select the pinned release workspace. In Terra's Data tab, select the required
   rows of `phenotype` and `participant`, then **Export → Download as TSV**.
   Export family context and molecular linkage tables needed for the selected
   use cases as well. A table export contains metadata and file references; it
   does not download the referenced VCFs.
   [Terra table export instructions](https://support.terra.bio/hc/en-us/articles/4417345161627-How-to-modify-and-edit-data-tables).
2. Preserve the untouched exports. Record workspace namespace/name, release,
   consent group, table names, export time, selection criteria, row counts,
   declared model version, and SHA-256 checksums in a protected acquisition
   manifest. Export from a stable release snapshot and verify cross-table
   references before publishing the bundle; several separate downloads are not
   an atomic multi-table snapshot.
3. Apply a versioned conversion for Terra entity-ID headers and reference-valued
   cells into the importer's canonical TSV contract. Record input/output hashes
   and conversion version; preserve original identifiers and values. Do not
   assume every exported cell is already in deposition-TSV format.
4. Build the observation mapping from release metadata linking participants to
   specimens/experiments and callset sample identifiers, then verify it against
   the VCF header. Record the actual tables/columns traversed for that release.
   No join based only on a filename or similar-looking ID is acceptable. Reject
   ambiguity; retain participants without genetic observations as clinical-only
   records and observations without clinical mappings as unmapped.

#### VCF discovery, download, and preparation

Use the release workspace's molecular-file references and release-specific
methods/errata to select the callset. GREGoR documents how to find file paths
through the Data tab and provides a programmatic notebook; do not hard-code a
bucket path inferred from a workspace name.
[Finding GREGoR key files](https://gregorconsortium.org/research-home/key-files).

For the first KING validation cohort, prefer one consistently processed joint
SNV callset with explicit sample genotypes, rather than combining heterogeneous
single-sample files. Public sites-only VCFs lack the required sample GT evidence
and cannot replace the controlled-access callset. Download each selected
`.vcf.gz` and its matching `.tbi` or `.csi`, recording object URI, generation,
size, available provider checksum, and local SHA-256. Check header sample sets,
reference declarations, and index usability before annotation.

The following is a template for an authenticated operator, not a runnable GREGoR
URL. Use manifest-resolved, generation-pinned objects and the authorized billing
project; repeat for the actual index object rather than assuming its suffix:

```bash
gcloud storage cp --billing-project=AUTHORIZED_PROJECT \
  'gs://AUTHORIZED_BUCKET/RELEASE/FILE.vcf.gz#GENERATION' \
  /approved/staging/FILE.vcf.gz
```

The CLI supports copying cloud objects and requester-pays billing options.
[Google Cloud storage copy reference](https://docs.cloud.google.com/sdk/gcloud/reference/storage/cp).

Stage transfers under temporary names/directories and promote them only after
verification. Reuse a completed local artifact only when its manifest and hashes
match; failed authentication, missing objects, interrupted transfers, or checksum
mismatches must not publish a usable bundle. Keep credentials and signed URLs
out of manifests and logs. Retain protected object locators only in the authorized
storage boundary. Do not download an entire bucket to obtain a small cohort.

Where practical, subset samples and regions inside the approved cloud environment
before localizing derived files, with the exact sample list, interval coordinates,
tool versions, and parameters recorded. Keep genome-wide panel coverage for
biological KING validation; a tiny locus slice tests plumbing only. Preserve
explicit reference and missing genotypes. Never interpret a missing variant-only
record as a reference call or feed unprocessed gVCF blocks into the KING loader.

Normalize/annotate the selected VCF upstream with pinned reference data and VRS
software, then prepare the exact panel and provenance headers required by the
[KING guide](how-to-king-robust.md). Do not invent reference calls or stamp headers
without verifying the preparation they describe. Record source and derived
hashes, annotation/reference versions, panel hash, sample renaming map, and QC
policy. Build the genetic index first, verify its artifact hash, then bind the
phenotype snapshot to it. Phenotype-only updates can reuse the verified genetic
artifact. Measure download, subsetting, annotation, and indexing separately from
query latency; network transfer is not a matcher benchmark.

### Atomic import and reproducibility

Add a cohesive `phenotypes.py` module for schema, import, and batched reads, with
CLI wiring in `cli.py`. Split it only if implementation size warrants it.
The illustrative interface is:

```text
vrs-matcher import-phenotypes --db cohort.db --phenotype-db clinical.db --manifest bundle.json
```

The manifest names participant, phenotype, observation-mapping, and optional
reconciliation files for each source. Support GREGoR deposition TSVs without
`phenotype_id` and table exports containing that field. Require participant
records for assertion references, even when optional family context is absent.

1. Pin the source release schema and its digest. Initially test both the story's
   reviewed 1.12 contract and a documented 1.11-to-canonical adapter for R05;
   retain source and target versions rather than relabeling the input. Validate
   TSV headers, required values, enumerations, references,
   and ontology/term namespace compatibility against that pinned contract.
   Preserve supported optional fields, modifiers, details, and source IDs.
   Reject unsupported model versions explicitly. Ontology membership checks
   require a pinned supplied vocabulary; never imply that prefix validation
   proves a term exists. No network resolution occurs during import or query.
2. Hash input bytes and a canonical manifest, excluding paths and import time
   from content identity. Derive `snapshot_id` from that manifest. Generate local
   assertion IDs by hashing source namespace, source-file digest, and row number;
   preserve the supplied phenotype ID separately. This keeps conflicting and
   repeated rows addressable. IDs are stable for the same input bytes, not for
   reordered files. Do not deduplicate using participant/term alone.
3. Stream validation and inserts in bounded batches under one transaction.
   Recheck input hashes before commit. Any invalid row, changed input, failed FK,
   ambiguous mapping, or interruption rolls back the complete new snapshot.
   An identical snapshot import is an idempotent no-op. Completed snapshots have
   no update/delete path in the initial importer.
4. Return counts of participants, assertions by presence, mapped/unmapped
   observations, repeated assertions, and conflicts. Errors include source and
   row identifiers without dumping clinical free text into logs.

Corrections create a new full snapshot. A snapshot-diff operation compares
scoped participant, ontology, term, presence, and normalized optional content,
not local assertion IDs alone; row reordering is not a biological change.
Keep missing optional values distinct from explicit values. Retain old snapshots
and original source artifacts so earlier reports remain reproducible. Import
time is provenance, not an inferred clinical assessment date.

### Query and output behavior

Add opt-in `--phenotype-db` and `--phenotype-snapshot` options to pair and
one-versus-all commands; require both together. Without them, preserve current
output and execution paths. With them:

1. Run the existing matcher, including its normal ranking and top-N behavior.
2. Collect distinct observation IDs in the returned results, including the query
   and unscorable results. Retrieve mappings and assertions with batched joins,
   deduplicating participant lookup across repeated observations.
3. Return a versioned envelope with the unchanged genetic payload, snapshot and
   genetic provenance, observation-to-person mappings, participant assertions,
   linkage status, and completeness/conflict summaries. Store each participant's
   assertions once in the envelope, referenced by results. Text output presents
   separate genetic and phenotype sections.

Use parameterized SQL with bounded batches (default 500 IDs, respecting the
connection's variable limit). No phenotype SELECT executes inside a scoring
loop. Do not parse the entire sidecar or cache it globally. Keep deterministic
ordering for assertions, source records, and participant summaries.

Presence is an assertion, not a forced single participant state. For a term,
report the set of recorded statuses and a conflict flag. An explicit `Unknown`
row differs from no assertion; no assertion differs from an unmapped observation.
A reconciled participant can have both positive and negative source assertions.
Count unique person keys, not observations or assertion rows. Report scoped
family counts and missing family information without assuming unknown families
are independent or that cross-source family IDs have been reconciled.

For all-pairs exports, write the genetic pair stream plus a separate participant
annotation stream and a manifest connecting them. Retrieve annotations in bounded
batches once per selected participant set, rather than repeating them for every
pair. Do not materialize all pairs to enrich them. The existing KING unscorable
list and existing sample-ID lists remain known memory costs; this feature must
not claim to bound those pre-existing structures.

### Explicit phenotype cohorts

Provide exact `(ontology, term_id, presence)` cohort selection, initially one
term predicate per request. A positive cohort includes people with at least one
explicit `Present` assertion, flags conflicting assertions, and offers an
explicit conflict-exclusion option. Report unknown-only, absent-only, conflicting,
and unrecorded groups with defined counts. No ontology expansion or semantic
ranking is implied. Repeated equivalent assertions never increase a score.

Resolve the cohort to observation IDs before scoring; intersect it with the
registered sample set. An explicit filter changes eligible peers, not marker
selection or the genetic equations. Apply top-N **after** restricting peers;
filtering an already truncated result is incorrect. Keep the named query
observation available even if it does not meet the peer predicate. Empty cohorts
return an empty result with cohort provenance.

For built-in matchers, use a small scoped `PluginContext` adapter whose
`list_samples()` returns sorted eligible IDs and whose other methods delegate
unchanged. Both built-ins already discover peers through this method, preserving
KING batching and heap selection. An all-pairs cohort contains only eligible
observations. Do not offer cohort restriction to third-party plugins until they
explicitly support this discovery contract; fail clearly rather than promising
a filter that a plugin can bypass. Post-score enrichment needs no protocol change.

Record predicate, conflict policy, snapshot, selected-ID digest and count, query
exception, and algorithm parameters. Never interpret phenotype selection as
independent evidence validating a genetic match.

## Performance contract

The sidecar separates clinical import I/O, storage growth, and writer locks from
the genotype database. It also makes phenotype updates independent of VCF parsing.
Its costs are snapshot duplication, explicit artifact management, and extra reads
only for requested enrichment. Simultaneous workloads can still contend for disk
and CPU; a separate file does not remove physical resource contention.

For R returned observations, P distinct linked people, and A retrieved assertions,
enrichment should use batched indexed lookups and O(R + P + A) output memory, not
scan the cohort's genotype rows or retain a pairwise phenotype matrix. Large
unlimited reports use streams. A top-N genetic result may still have many
unscorable peers and large phenotype payloads; measure and report those costs.

Use `EXPLAIN QUERY PLAN` and SQL tracing on representative populated sidecars to
verify indexed cohort and participant lookups and absence of per-pair phenotype
queries. Avoid brittle assertions on SQLite's exact plan text. Measure index
build overhead and storage before adding any additional indexes. Defer packed
phenotypes, denormalization, multiprocessing, and caching until profiling shows
a bottleneck.

## Test strategy

Use synthetic clinical records only in committed fixtures. Extend the existing
pytest suite and benchmark tooling; no external GREGoR access is required in CI.

| Test layer | Required cases and assertions |
|---|---|
| Import and constraints | TSV with/without source phenotype IDs; all optional fields preserved; missing required values, invalid enums, bad references and unsupported schemas rejected; same local IDs in separate namespaces remain separate. |
| Mapping/reconciliation | Two observations count as one person; ambiguous observation mapping fails; explicit cross-source reconciliation exposes both sources; unreconciled IDs never merge; unmapped samples remain visible. |
| Assertion semantics | Present/Absent/Unknown/no row/unmapped remain distinct; conflicts survive; repeated rows do not inflate participant counts; namespaces and broader/narrower terms are not treated as equivalent. |
| Atomicity/replay | Inject a failure after an insert batch and before commit; no partial snapshot is visible and older snapshots remain readable. Reimport is idempotent; file mutation fails; reordered rows do not appear as clinical changes in a snapshot diff. |
| Genetic regression | Pair, one-versus-all, ties, candidate-marker restrictions, empty calls, and unscorable cases retain identical genetic fields and order with enrichment enabled. Genotype-only JSON/text remains compatible. |
| Cohort correctness | Compare scoped built-in results with an exhaustive pairwise oracle over exactly the eligible peers, then sort/truncate; include a globally top-ranked excluded peer, empty cohorts, query outside cohort, and conflicts. Unsupported custom-plugin filtering fails explicitly. |
| Reanalysis/provenance | Changing only phenotypes leaves genetic scores unchanged; old snapshots reproduce old reports; wrong/missing artifact binding fails validation; mapping corrections create new snapshots. |
| Performance structure | Trace SQL to prove zero sidecar opens/queries without opt-in; lookup count grows with batches, not pairs; all-pairs export consumes iterators and emits participant annotations separately. |
| Scientific output | Reports never label an allele causal, infer identity from phenotype agreement, infer absence from missing data, or present carrier overlap as complete segregation evidence. |

A curated, authorized analyst review follows automated tests. Include relatives,
replicate observations, sparse records, discordant assertions, and different term
granularity. Review linkage accuracy, usefulness for case review, and missed
cohort members; do not treat clinical agreement as genotype ground truth.

### Fixture ownership, versioning, and refresh

Maintain three distinct fixture tiers; their locations and commands below are
proposed implementation work, not existing artifacts:

| Tier | Storage and contents | Execution and maintenance |
|---|---|---|
| Offline regression | Commit small, wholly synthetic TSVs, raw/annotated VCFs, panels, mapping files, manifests, and expected semantic results under `tests/data/gregor/`. | Default pytest; no credentials, network, cloud CLI, or SeqRepo download. Cover both supported schema versions and Terra/deposition formats. |
| Synthetic scale | Commit a deterministic generator and configuration/seed; generate large TSVs/VCFs in temporary benchmark storage. | Dedicated benchmark job; record generator revision, versions, hashes and sizes. Do not commit generated databases or large cohorts. |
| Authorized integration | Keep a minimal real-data subset, manifests, expected linkage review, and derived artifacts in access-controlled project storage outside the repository. | Explicit opt-in job in an approved environment; no public CI downloads or artifact uploads. Access and retention follow the source authorization. |

The offline fixture README must state that people, families, clinical text, and
identifiers are invented, document each scenario, and distinguish synthetic VRS
IDs used for loader tests from IDs actually computed for a small artificial
reference. Include a small reference/annotation fixture or an optional pinned
annotation integration check so fabricated identifiers are not mistaken for
biological validation. Reuse existing KING oracle cases where applicable.

Maintain at least two phenotype snapshots over the same synthetic genetic input,
plus conflicting sources, reconciliation, repeated observations, missing mappings,
all presence states, reference/missing GTs, and malformed-input cases. Expected
results should be hand-reviewed counts, mappings, and genetic oracle values;
never regenerate them automatically from the implementation under test. Compare
semantic content separately from byte hashes when compression metadata or tool
versions affect derived bytes.

Provide one documented fixture regeneration command with fixed seeds and tool
versions. A refresh PR must include its reason, changed input hashes, schema or
conversion changes, expected-result differences, and review by a maintainer
familiar with the biological scenarios. The implementer owns initial fixtures;
subsequent importer/schema changes must update them in the same PR. Keep older
supported-version fixtures to prevent regressions. A new GREGoR release triggers
compatibility review, not automatic replacement or a background download.

For restricted integration runs, require an explicit enable flag and a manifest
path (for example `VRS_MATCHER_GREGOR_INTEGRATION=1` and
`VRS_MATCHER_GREGOR_MANIFEST`). Without opt-in, report skipped tests. Once enabled,
missing files, unavailable authorization, and hash mismatches fail the run rather
than silently skipping validation. The fetch step is separate from pytest;
provide regeneration/subsetting recipes in protected storage and retain private
truth labels and participant mappings there. Renaming real participants does not
make their genotypes or clinical records suitable for committing as synthetic
fixtures. If access expires, follow the required retention/deletion policy and
record the resulting reproduction limitation; never copy the fixture into a
public location to preserve a test.

Add offline acquisition tests using local files and stubbed command responses:
TSV header conversion, sample linkage, wrong release/model, missing VCF index,
truncated copy, checksum mismatch, interrupted staging, and idempotent reuse.
Assert that no incomplete bundle reaches ingestion. CI should check the fixture
manifest checksums and absence of credentials, signed URLs, or real-data exports
in tracked fixture files; review remains necessary for clinical content.

### Benchmark and release gates

Compare the parent revision with the implementation on identical genetic inputs
and hardware. Reuse the [KING benchmark phases](adr-king-performance.md), adding
separately measured phenotype import, cohort lookup, enrichment, serialization,
and snapshot-diff phases. Measure identity queries too.

Use deterministic fixtures spanning 1,000, 10,000, and 100,000 participants with
10 and 100 assertions per participant, plus skewed high-annotation participants,
replicate mappings, conflicts, and sparse coverage. Vary the size of the phenotype
sidecar independently of a fixed genetic cohort to expose accidental full scans.
Use representative genetic marker counts for query runs. Bound all-pairs trials
to practical subcohorts; do not generate a 100,000-person pair matrix merely to
test phenotype scaling.

Record input hashes, versions, hardware, batch size, output size, wall/CPU time,
peak RSS, database size, query counts, and plans. Report first-run measurements
separately, then median and range over at least five warm repetitions. State
filesystem-cache conditions; process restart alone is not a cold-disk test.

Proposed gates, to be measured rather than asserted as current performance:

- Exact genetic result parity and zero phenotype I/O on genotype-only paths are
  mandatory deterministic checks.
- On a fixed reference runner, genotype-only median time and peak RSS must stay
  within 5% of baseline for workloads long enough to exceed timer noise. Repeat
  suspected regressions; investigate persistent failures before release. Timing
  gates belong to a controlled benchmark job, not ordinary shared-runner unit CI.
- Enrichment must not add genotype retrieval or alter scorer timings materially.
  At fixed returned participants, increasing unrelated sidecar data tenfold must
  not cause tenfold enrichment time or memory growth; inspect plans and traces
  for a full scan if it does. Report absolute latency and output volume rather
  than promising a universal enrichment percentage.
- Import and streaming-export transient memory must follow batch and record
  size, not total assertions or total pair count. Verify with multiple data
  sizes and unusually large detail fields; account separately for SQLite caches,
  existing matcher memory, and intentionally materialized output.

## Alternatives and consequences

Embedding clinical tables in the genetic database would simplify cross-table
FKs, but couples clinical writes and file growth to a large reusable genotype
artifact. The sidecar accepts explicit binding validation to avoid that coupling.

Joining phenotype rows into genotype scoring would multiply rows, risk changing
counts, and repeat clinical work per pair. Post-score enrichment avoids those
failure modes. A combined phenotype/genotype ranking plugin is deferred because
the story requires interpretable evidence, not an unvalidated composite score.

Full snapshots consume more disk than deltas but make reproduction, rollback,
and source conflicts straightforward. Delta storage and ontology reasoning can
be designed later if measured needs justify them.

## Implementation sequence and completion

1. Document and validate release acquisition, schema adapters, and fixture
   regeneration. Implement sidecar schema, pinned input contract, atomic importer,
   mapping/reconciliation validation, and snapshot diff with focused tests.
2. Add batched enrichment and versioned report output; verify genetic parity and
   unchanged default CLI/plugin behavior.
3. Add exact cohort selection, scoped built-in traversal, streaming exports, and
   the cohort oracle tests.
4. Run the scaling benchmarks and analyst review; document measured limits and
   usage examples. Keep status Proposed until the decision is adopted, and mark
   implementation completion separately from scientific validation.

Completion requires the user story's acceptance scenarios, the deterministic
regression checks, reproducible performance evidence, and documentation of
artifact binding and authorized clinical-data handling. No existing index needs
migration and no genetic scoring algorithm changes as part of this decision.
