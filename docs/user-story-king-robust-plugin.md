# User story: Manichaikul (KING) robust relationship estimation plugin

Status: Core implemented · 2026-09-30. Acceptance A–D are tested; E has an
opt-in benchmark/report runner and a recorded synthetic smoke run, but biological
validation is pending a genome-wide panel with verified pedigree/replicate labels.
See the [usage guide](how-to-king-robust.md) and [validation report](king-validation.md).
The specification below retains the intended biological acceptance requirements;
the complete definition of done is not yet satisfied by a synthetic fixture.

## User story and value

As a bioinformaticist investigating sample identity and familial relationships, I want to estimate pairwise kinship from a VRS-indexed SNP panel and rank potential relatives, so that I can distinguish duplicate-like observations from other close relationships and inspect the genotype evidence behind each estimate.

The plugin should reuse the repository's indexed sample workflow while providing a recognized relationship estimator. Its validation should establish numerical correctness, biological usefulness, and the cost of indexing versus repeated analysis. It must not reinterpret the existing Jaccard or weighted-concordance scores as kinship.

This story extends the [VCFtools benchmark ADR](adr-benchmark-vcftools.md) with a separate relatedness track. The ADR's genotype-discordance baseline remains appropriate for its original task; this track uses `vcftools --relatedness2`.

## Scope and estimator contract

Implement a plugin named `king-robust` for human autosomal, diploid, biallelic SNVs on one declared reference and a fixed, versioned panel. Input VCFs must contain valid, consistently generated VRS ALT identifiers. Exclude multiallelic records, indels, symbolic alleles, sex chromosomes, mitochondrial variants, and non-diploid calls in this first version. Do not silently split multiallelic records into independent markers.

For each pair, calculate counts only over panel loci with a passing, complete diploid genotype in **both** observations:

| Symbol | Meaning |
|---|---|
| M | Number of jointly called eligible loci |
| H_i, H_j | Heterozygous calls in each individual on those loci |
| HH | Loci where both individuals are heterozygous |
| O | Opposite homozygotes: `0/0` versus `1/1`, in either direction; IBS0 count |

Use the between-family estimator for ranking unknown relationships:

```text
kinship = 0.5 - (H_i + H_j - 2*HH + 4*O) / (4*min(H_i, H_j))
```

Also return the within-family estimate as an explicitly named diagnostic:

```text
kinship_within_family = (HH - 2*O) / (H_i + H_j)
```

These correspond to equations 11 and 9, respectively, in Manichaikul et al. The paper uses them for between-family and within-family checking. This first plugin applies the between-family estimator by default because the index has no pedigree assignments; it does not implement pedigree-aware switching or the full KING application. [Original paper](https://www.cog-genomics.org/static/pdf/Manichaikuletal2010.pdf)

Preserve negative estimates. If a formula's denominator is zero, return a null value with an explicit reason; never substitute zero kinship. For M >0, also return `ibs0_fraction = O/M`; otherwise it is null. A valid numeric estimate on a tiny panel is not a validated relationship call.

VCFtools' inspected implementation uses the within-family expression, but accumulates heterozygote totals per individual rather than restricting them to jointly called loci. Therefore, require direct parity for the within-family diagnostic on fully called panels, and test missing-data differences explicitly. Pin the executable and source revision used; do not label the default between-family score as identical to VCFtools. [VCFtools implementation](https://github.com/vcftools/vcftools/blob/master/src/cpp/variant_file_output.cpp#L4353)

## Required changes to indexing

### Preserve called reference genotypes

The current [loader](../src/vrs_matcher/loader.py) emits carried ALT alleles only, skips reference and fully missing genotypes, and classifies a partial call such as `./1` as HET. The current [database](../src/vrs_matcher/db.py) cannot reconstruct callability from absent allele rows. These behaviors cannot supply unbiased KING counts on general inputs.

Add an opt-in genotype-indexing path to ingestion. Keep the carried-allele table and existing identity behavior intact. A minimal additional SQLite representation is:

| Record | Required content and constraints |
|---|---|
| Panel metadata | Panel identifier/hash, reference identity, annotation version, and complete eligible marker list. One panel per genotype index in this version. |
| Panel marker | Unique ALT VRS ID, reference sequence identity, position, REF, ALT; enforce one biallelic marker per locus and consistent allele orientation. |
| Indexed observation | Unique sample/observation ID, source provenance and input hash, QC policy, completed genotype-index capability/version, and ingestion counts. Register observations even when they have no passing calls. |
| Genotype call | Unique `(sample_id, vrs_id)` with ALT dosage constrained to 0, 1, or 2, referring to a panel marker and indexed observation. Include **dosage 0** for passing `0/0` calls. |

Only passing, fully called genotypes need call rows. The presence of a row establishes callability; absence means unknown or excluded, **never homozygous reference**. Thus explicit rows for every missing cell are unnecessary. Retain ingestion summaries by exclusion reason so absent annotation, failed QC, missing GT, and unsupported records remain auditable. No imputation is performed.

### Ingestion rules

1. Select the declared panel before writing genotype rows. Resolve contig aliases against the declared reference; chromosome labels alone do not establish sequence identity. Reject inconsistent marker identity, REF/ALT orientation, or ambiguous duplicate panel loci. Do not infer reference calls for panel sites missing from a variant-only VCF or expand gVCF reference blocks in this version.
2. Require exactly one valid ALT VRS identifier for each eligible biallelic panel record. Reject missing/placeholder annotations, wrong annotation cardinality, and marker mismatches. VRS generation and reference normalization remain upstream responsibilities; record their provenance.
3. Apply PASS/unfiltered record policy and configurable GQ/DP thresholds, initially matching existing defaults GQ ≥20 and DP ≥0 when values exist. Missing QC values remain allowed and counted. Use the same filters for dosage 0, 1, and 2. For the KING genotype table, partial GT, fully missing GT, invalid allele indices, and non-diploid GT produce no call row. Phase and allele order do not change dosage.
4. Give independently measured or processed observations distinct sample IDs before loading. `source_dataset` does not namespace the current primary key. Reject reuse of an observation ID in this initial genotype-indexing path; require an explicitly rebuilt index to replace an observation. Never silently merge or overwrite conflicting calls.
5. Reject duplicate panel records within an observation rather than allowing order-dependent replacement. Publish calls and the completed observation marker atomically; a failed/interrupted load must not leave a usable partial genotype index. Batch inserts may occur within that transaction.
6. Return counts for indexed observations, retained markers, dosage 0/1/2 calls, and exclusions by reason. Store the panel and QC manifest with the database. Use a single consistent QC policy per genotype index; reject incompatible append settings.

Existing databases remain usable by `identity`. A KING request against a database or sample lacking the completed genotype-index capability must fail with a clear re-ingestion instruction. Adding empty tables alone must not mark old data as KING-ready. The original VCFs are required to populate the new data; reference and missing calls cannot be recovered from the old index.

## Plugin, API, and CLI behavior

- Extend [PluginContext](../src/vrs_matcher/plugins.py) with read-only access to panel/capability metadata and called genotypes, including dosage 0. Compute on the pairwise intersection of called markers, not just shared carried ALT IDs.
- Allow `candidate_vrs_ids` to select markers from the declared panel; recalculate all counts on the restricted panel and reject unknown marker IDs. An empty restriction yields no comparable loci, not a perfect match.
- Add a distinct kinship result type and corresponding CLI rendering. Return sample IDs, estimator name/version, panel ID, M, H_i, H_j, HH, O, both estimates, IBS0 fraction, and status/reason. Preserve full precision in machine-readable output. Do not place kinship in fields labeled Jaccard or weighted concordance.
- Adapt the plugin return contract and dispatch to accept the new result type while preserving existing identity results and local/entry-point plugin behavior. If a breaking protocol change is unavoidable, version it explicitly and report unsupported plugins clearly.
- Rank finite between-family estimates in descending order, breaking ties deterministically by sample ID. Exclude the query itself from one-versus-all results. Report unscorable pairs separately, outside top-N. Pairwise queries return their explicit unscorable result. `shared-variants` must reject kinship results clearly unless actual shared-allele output is separately supported.

Illustrative CLI, to be implemented as part of this story:

```bash
vrs-matcher load-samples cohort.vrs.vcf.gz --db cohort.db \
  --index-genotypes --panel king-panel.tsv
vrs-matcher match-samples SAMPLE_A SAMPLE_B --db cohort.db --algorithm king-robust
vrs-matcher match-sample SAMPLE_A --db cohort.db --algorithm king-robust --top 10
```

Return estimates and evidence counts in this version; automated relationship labels and pedigree reconstruction are outside scope. Kinship does not distinguish parent–offspring from full siblings by itself. Retaining IBS0 supports subsequent interpretation. Population admixture can also bias estimates, so report population-stratified validation rather than promising universal robustness. [KING estimator documentation](https://www.cog-genomics.org/plink/2.0/distance)

## Acceptance tests

Implement deterministic tests in the existing pytest suite. Keep external-tool, public-data, and performance runs opt-in. Each test must exercise the new implementation; documentation of expected behavior alone does not satisfy acceptance.

### A. Indexing and data integrity

| ID | Given / when | Required outcome |
|---|---|---|
| A1 | Passing `0/0`, `0/1`, `1/1`, `./.`, and `./1` calls are ingested. | New table contains dosage 0, 1, and 2 only; missing and partial calls are excluded and counted. The original identity index retains its established behavior. |
| A2 | One sample has `1/1`; its peer is `0/0`, missing, filtered, or absent at the marker in four separate cases. | Only explicit passing `0/0` contributes an opposite-homozygote count. Other cases do not contribute to M or O. |
| A3 | GQ values 19, 20, and missing; DP below and at a configured positive threshold; PASS, unfiltered, and failed FILTER. | Boundary inclusion matches policy for all dosages, including reference calls; exclusion summaries reconcile with input counts. |
| A4 | Inputs include indels, multiallelic sites, non-autosomal loci, haploid/polyploid calls, invalid allele indices, and reversed phased GT order. | Unsupported records/calls are excluded and counted; invalid or missing alleles never become dosage 0; valid diploid phase/order changes preserve counts. |
| A5 | Missing IDs, placeholder IDs, inconsistent marker orientation/reference, duplicate records, reused observation IDs, or incompatible append policy. | Actionable error; no silently overwritten calls or newly published partial observations. Inject a mid-load failure and verify rollback. |
| A6 | An old allele-only database, or a mixed database with a sample not genotype-indexed, is queried. | KING refuses the unsupported sample/database with a re-ingestion message; identity still works. A completed all-missing observation instead returns an unscorable result. |
| A7 | Round-trip a valid genotype index, close it, and reopen it. | Metadata, reference calls, dosages, summaries, and scores are unchanged. |

### B. Hand-calculated statistical oracle

Use this eight-marker fixture, with dosage 0=`0/0`, 1=`0/1`, and 2=`1/1`:

| Markers | Sample i dosage | Sample j dosage |
|---|---|---|
| 1–3 | 1 | 1 |
| 4 | 1 | 0 |
| 5–7 | 0 | 1 |
| 8 | 0 | 2 |

Expected counts: M=8, H_i=4, H_j=6, HH=3, O=1. Expected between-family kinship **0.0**, within-family estimate **0.1**, and IBS0 fraction **0.125**. Assert integer counts exactly and floating-point values within `1e-12`.

Additional required cases:

| ID | Case | Required outcome |
|---|---|---|
| B1 | Identical observations with at least one heterozygous marker, including a self-comparison. | Both estimates 0.5 and O=0. |
| B2 | Two markers: `(1,1)` and `(0,2)`. | Both estimates −0.5; no clipping to zero. |
| B3 | No common calls, all reference calls, or one sample with zero heterozygotes. | M and other counts remain truthful; between-family estimate is null when its denominator vanishes. Within-family estimate is independently null only when its denominator vanishes. No NaN/Infinity output. |
| B4 | One shared `(1,1)` marker plus one `(1,missing)` marker. | M=1, H_i=H_j=HH=1, O=0; both estimates remain 0.5. A heterozygote at a non-jointly-called marker must not alter the denominator. |
| B5 | Swap sample order, reorder records, change phase, or consistently reverse REF/ALT coding in both samples in a separately valid panel. | Estimates and symmetric counts are invariant; H_i and H_j swap when sample order swaps. |
| B6 | Restrict the eight-marker fixture to markers 1–3, then to an empty set. | First result has M=3 and both estimates 0.5; second has M=0 and null estimates. |

### C. Plugin integration and compatibility

Plugin discovery lists `king-robust`. Pair and one-versus-all API/CLI paths expose correctly named scores and counts. Compare every ranked result to a direct pair calculation; verify deterministic ties, top-N behavior, self exclusion, and separate unscorable reporting. Confirm candidate restrictions include explicit reference calls. Existing identity, plugin, loader, and CLI tests continue to pass. No changed Jaccard semantics or reference rows added to the carried-allele API are permitted.

### D. VCFtools numerical comparison

Pin and record VCFtools version/build and run `vcftools --gzvcf panel.vcf.gz --relatedness2 --out baseline` on the same prefiltered, fully called biallelic SNP panel used for indexing. Materialize QC decisions before either tool executes. Compare every pair's `N_AaAa`, `N_AAaa`, `N1_Aa`, and `N2_Aa` exactly against HH, O, H_i, and H_j. The plugin's **within-family** diagnostic must match `RELATEDNESS_PHI` within half the last printed decimal unit, allowing a small floating-point rounding tolerance.

Do not require default between-family parity: the eight-marker oracle intentionally gives 0.0 versus VCFtools' 0.1. Include B4 as a missingness diagnostic. In the inspected VCFtools implementation, its per-individual heterozygote totals would produce 1/3 rather than the plugin's 0.5. Verify and report the pinned executable's actual behavior; investigate deviations instead of changing the plugin's pairwise callability rule to match them. [Source of comparator semantics](https://github.com/vcftools/vcftools/blob/master/src/cpp/variant_file_output.cpp#L4353)

### E. Biological validation and performance report

Use a frozen genome-wide SNP panel with documented duplicate/replicate, parent–offspring, sibling, and unrelated labels. Keep donor families separated between any calibration and evaluation partitions. Exact duplicates are positive controls, not substitutes for independently measured replicates. The existing chr22 slice is only a smoke test.

Deliver pair-level estimates/counts and relationship-stratified summaries, duplicate retrieval accuracy with ties reported, and an error analysis stratified by missingness and ancestry. Report sample/pair counts and uncertainty; compare the two named estimators separately. Statistical accuracy requires external labels, not agreement with VCFtools. No universal clinical classification threshold or minimum sensitivity is claimed by this story; such deployment criteria require a separately calibrated decision rule.

Follow the [benchmark ADR's measurement protocol](adr-benchmark-vcftools.md): record preparation, VRS annotation, genotype/allele indexing, single-query and all-pairs wall time, peak memory, storage, versions, hashes, hardware, and repeated-run variability. Compare identical pair sets; VCFtools' all-pairs output must not be presented as a single-pair cost comparison. Include reference-call storage overhead and failures at larger sizes. Completion requires reproducible measurements and explained discrepancies, not an assumed speedup or statistical improvement.

## Definition of done

The feature is complete when A–D pass, the opt-in E report is reproducible from pinned inputs, the indexing and estimator contracts are documented in usage help, and existing identity/plugin behavior remains compatible. Deliver one implementation and the smallest tests needed to exercise these cases. Defer bitset optimization, alternative storage engines, cross-assembly matching, automatic pedigree inference, and additional variation classes until a measured need justifies them.
