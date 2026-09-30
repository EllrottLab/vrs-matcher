# ADR: Separate KING benchmark workloads and optimize matching incrementally

## Status

Accepted · 2026-09-30. The benchmark reports index construction, one-versus-all
query, indexed all-pairs matching, and direct VCFtools execution as distinct
workloads. Batched genotype retrieval, pair-result streaming, and bounded top-N
selection are implemented. A deterministic synthetic benchmark mode now reports
index construction, SQLite call retrieval, and scoring/iteration separately.
Packed genotype representations and native or parallel scoring remain
conditional follow-up work.

## Context

The KING benchmark combines work with different lifecycles: building the SQLite
index, answering a single query against an existing index, matching all pairs
from an existing index, and running VCFtools directly on a prepared VCF. A single
performance percentage hides index cost and can imply equivalence where the
workloads differ.

The current all-pairs plugin repeatedly retrieves each sample's genotype-call
dictionary through `PluginContext` for every pair. For a cohort of N observations,
many of the same calls are read and converted repeatedly. Pair scoring then
allocates a set intersection and traverses it in Python. These are candidate
bottlenecks, but the existing tiny smoke fixture is dominated by process startup
and cannot establish their importance at realistic scale.

## Decision

### Benchmark workloads

Keep the isolated timing rows and report workload comparisons separately:

| Workload | Includes | Use |
|---|---|---|
| Index build | VRS-Matcher genotype and allele SQLite ingestion | Quantify initial indexing cost |
| One-versus-all query | Query against an already-built index | Quantify repeated single-query latency |
| Indexed all-pairs | Pair scoring against an already-built index | Isolate VRS-Matcher matching after ingestion |
| Index + all-pairs | Index build and all-pairs work in paired repetitions | Show VRS-Matcher one-shot pipeline cost |
| VCFtools direct | Read the prepared VCF and run `--relatedness2` | External baseline for the current direct-file workflow |

Report first run separately from subsequent-run median/range. Calculate medians
for composed workloads from per-repeat sums, pairing each VRS-Matcher repeat
with that repeat's VCFtools run. The relative percentage remains
`100 × VCFtools median time / VRS-Matcher median time`, with VCFtools at 100%.
Label it as a relative timing ratio, not a universal speed claim.

Neither comparison is fully like-for-like: indexed all-pairs excludes VCF
reading for VRS-Matcher but not for VCFtools; index + all-pairs includes
VRS-Matcher indexing while VCFtools reads an already-prepared VCF. The two
programs also emit different pair sets. Preserve these caveats in every report.
Do not invent a VCFtools single-query comparison; the tool does not provide that
workload in this runner.

### Matcher architecture improvements

Preserve SQLite as the durable index and preserve the existing scoring equations
and result contract. Before changing representations:

1. **Implemented:** `PluginContext.get_called_genotypes_many` fetches calls in
   SQL batches. All-pairs matching processes sample tiles of 128, loads no more
   than two call tiles at once, and yields scores as an iterator. The benchmark
   writes the pair TSV incrementally. Pair scoring remains a pure in-memory
   function.
2. **Implemented:** one-versus-all matching loads peer calls in batches and
   retains only top-N scorable results in a bounded heap. Ordering is
   deterministic. Unscorable results remain fully reported as required by the
   existing result contract.
3. **Deferred:** profile representative cohorts for the memory/time tradeoff
   of tile size and for the cost of rereading calls across sample tiles. The
   default tile size bounds transient genotype maps, but all unscorable results
   are still retained in memory to preserve current behavior.
4. **Implemented instrumentation:** `scripts/benchmark_king.py` can generate
   deterministic synthetic cohorts (`--synthetic-samples`, `--synthetic-markers`,
   `--seed`) and records internal index-build, tiled retrieval, and scoring plus
   pair-iteration durations separately. It also records sample/marker/pair
   counts and input hashes. Profiling retrieval and scoring uses the same
   bounded tiled iterator as the production all-pairs path.
5. Only if profiles still show Python dictionary/set traversal as a material
   cost, prototype packed dosage arrays or per-dosage bitsets for the called
   genotype matrix. Keep the SQLite tables as the canonical persisted format
   initially; build and cache the compact representation at query time or in an
   explicitly versioned derived index. Validate all counts and kinship values
   against the existing implementation before considering a storage migration.
6. Consider native/vectorized scoring or parallel execution only after the
   batch and representation experiments identify the remaining bottleneck.

## Consequences

- Reports distinguish index amortization from matching latency and avoid
  presenting a single percentage without a workload definition.
- Batch retrieval removes repeated per-pair SQL lookups and Python object
  construction; tiled iteration bounds transient call-map memory and streaming
  bounds retained pair-output memory.
- Compact representations may improve dense-panel scoring but add complexity,
  conversion cost, and memory overhead. No such representation should be
  adopted based only on the eight-marker smoke fixture.
- Existing SQLite databases and result semantics remain compatible through the
  first optimization step.

## Validation and rollout

An in-process synthetic comparison used 250 observations, eight markers, and
31,125 unordered pairs. Three timed runs per path produced median all-pairs
times of 1.892 s for the previous `match_pair` API loop and 0.154 s for bounded
batched iteration (12.3x faster, 91.9% lower elapsed time). Every score/count
matched exactly. The tiny eight-marker panel is not representative of
biological throughput, and this run did not measure peak RSS.

Benchmark fixed, reproducible cohorts at multiple sample and marker counts.
Record index time, one-query time, all-pairs time, total pipeline time, CPU,
peak RSS, and index size separately. Include cold process startup as the current
runner does, and explicitly state that filesystem cache is uncontrolled.
Compare exact counts and results against the current scorer, including missing
calls, candidate restrictions, unscorable pairs, ties, and top-N ordering.
Retain a compact smoke benchmark for report generation only; do not use it as
evidence of a cohort-scale speedup.

The synthetic mode isolates cost centers but does not model real ancestry,
missingness, annotation, or VRS identifier derivation. Use it to compare
implementation revisions on identical generated inputs, then validate scaling
on representative prepared data before making biological-workload performance
claims.

Initial phase-instrumentation run: deterministic synthetic cohort, 128
observations × 512 markers, seed 0, Python 3.13.5, three total runs (one first
run plus two subsequent runs). The subsequent medians were 0.639 s for index
construction, 0.111 s for tiled call retrieval, and 0.766 s for scoring plus
pair iteration over 8,128 unique pairs. The profile worker's total was about
0.877 s median, excluding startup. This is a repeatability and instrumentation
check, not a biological throughput result; the limited repeat count and
uncontrolled OS cache should be considered when interpreting it.

## Alternatives considered

- **Keep only one end-to-end percentage:** rejected because it conceals the
  index-build versus repeated-query tradeoff.
- **Compare indexed matching time directly as if equivalent to VCFtools:** not
  claimed; VCFtools still reads and parses the input VCF.
- **Replace SQLite with bitsets immediately:** deferred until representative
  profiling demonstrates that storage/retrieval or dictionary traversal is the
  dominant cost and memory tradeoffs are measured.
- **Add multiprocessing first:** deferred because it may increase memory and
  coordination overhead before repeated retrieval and Python-level traversal
  are addressed.
