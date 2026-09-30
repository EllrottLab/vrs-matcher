# KING implementation validation

Date: 2026-09-30. This report records implementation validation and a synthetic
smoke benchmark, not biological accuracy or a production performance claim.

## Acceptance status

| Area | Evidence | Status |
|---|---|---|
| A: indexing and integrity | Actual cyvcf2 VCF parsing, reference-call preservation, QC, exclusions, panel validation, rollback, append protection, and reopen tests | Passed |
| B: statistical oracle | Hand-calculated counts and equations; missing/partial calls, zero denominators, negative estimates, phase/orientation invariance, and candidate restrictions | Passed |
| C: integration | API/CLI, full-precision JSON, discovery, ranking/ties/unscorable peers, and existing identity/script-plugin tests | Passed |
| D: VCFtools comparison | External VCFtools 0.1.16, all pair directions including diagonals, fully called and missing-call fixtures | Passed |
| E: biological validation | No genome-wide panel with independently verified replicate/pedigree labels supplied | Pending; the runner is implemented and smoke-tested |
| E: performance | Six fresh-process repetitions on the synthetic two-observation, eight-marker fixture | Smoke only; scaling and end-to-end annotation costs remain unmeasured |

Validation commands:

```bash
uv run ruff check .
uv run ruff format --check .
uv run pytest --no-cov -q
VCFTOOLS=/path/to/vcftools uv run pytest \
  tests/integration/test_king_vcftools.py --run-integration --no-cov -q
```

Results: **134 passed, 3 opt-in tests skipped** in the ordinary suite;
**2 passed** in the explicitly enabled VCFtools test. The unrelated existing
1000 Genomes network/SeqRepo integration test was not run. Ruff lint and format
checks passed. No new runtime dependencies were added.

## Comparator and numerical results

VCFtools was built from the upstream `v0.1.16` source archive. The local
autotools launcher referenced an unavailable Perl executable, so the same C/C++
sources were compiled directly using the system compiler and zlib:

```bash
cc -O2 -c src/cpp/bgzf.c src/cpp/knetfile.c
c++ -O2 -std=c++11 '-DPACKAGE_VERSION="0.1.16"' \
  src/cpp/*.cpp bgzf.o knetfile.o -lz -o vcftools
./vcftools --version
```

The output was `VCFtools (0.1.16)`. No VCFtools source modifications were made.
The upstream archive is available at
[the v0.1.16 release source](https://github.com/vcftools/vcftools/archive/refs/tags/v0.1.16.tar.gz).

Archive SHA-256: `575c13073efe65cbff6e2ab99eef12fe04536f5dc1f98de6674c848ca83cf302`.

Built executable SHA-256: `925e84a70bb34da6f0dbf77016f9c9950f237c28ddecf6e1c0bfcaa9bd070641`.

| Case | Plugin between-family | Plugin within-family | VCFtools phi |
|---|---:|---:|---:|
| Eight-marker A/B oracle (M=8, H_i=4, H_j=6, HH=3, O=1) | 0.0 | 0.1 | 0.1 |
| One shared heterozygote plus one heterozygote missing in the peer | 0.5 | 0.5 | 0.333333 |

The second row is an expected semantic difference: VCFtools includes the
unpaired heterozygote in its individual total, while the plugin counts only
jointly called markers. The between-family and within-family formulas also
remain deliberately distinct. Detailed tests assert every underlying count,
with floating comparison tolerance derived from VCFtools' printed precision.

## Reproduce the smoke measurement

```bash
uv run python scripts/benchmark_king.py \
  --vcf tests/data/king/oracle.vcf --panel tests/data/king/panel.tsv \
  --vcftools /path/to/vcftools --output /tmp/king-smoke --smoke
```

The output directory must be new. Inputs contain synthetic identifiers and are
not suitable for biological inference. The complete output includes raw logs,
per-pair scores, comparator results, input/source/executable hashes, process CPU
and RSS measurements, and SQLite artifacts. Machine load and OS cache were not
controlled; first runs are not labeled cold-cache. Measurements include startup.
This final smoke run was performed after the test processes had completed.

Platform: `macOS-15.8-arm64-arm-64bit`; architecture: `arm64`; CPUs reported: 10.

Benchmark environment at time of measurement (before the Python 3.13 migration):
Python `3.12.1 (v3.12.1:2305ca5144, Dec  7 2023, 17:23:38) [Clang 13.0.0 (clang-1300.0.29.30)]`;
SQLite `3.43.1`; cyvcf2 `0.32.1`.

Unique pairs: 1.

Input preparation and VRS annotation are excluded; this is not end-to-end timing.
VCFtools computes both pair directions and diagonals; the plugin all-pairs stage
computes each unordered distinct pair once. Query is a separate workload.

Workloads are reported separately: `load` builds the genotype and allele
SQLite indexes; `query` matches one sample against the cohort; `all` matches
all unordered distinct pairs using an existing index; `vcftools` reads the
prepared VCF and computes `--relatedness2`. First-run time is repeat 0.
Subsequent median and range exclude repeat 0. Peak RSS is the maximum process
peak across all repeats for that stage, in MiB (1024² bytes). All times include
process startup.

| Stage | First run (s) | Subsequent median (s) | Subsequent range (s) | Max peak RSS (MiB) |
|---|---:|---:|---|---:|
| load | 0.168764 | 0.163348 | 0.158645–0.166269 | 37.80 |
| query | 0.162552 | 0.152377 | 0.149840–0.155873 | 36.30 |
| all | 0.153747 | 0.157489 | 0.150760–0.200823 | 36.23 |
| vcftools | 0.014874 | 0.006001 | 0.004983–0.007961 | 2.12 |

## Workload comparisons

Relative performance is `100 × VCFtools median / VRS-Matcher median`; VCFtools
is 100%. Each VRS-Matcher median is computed from paired subsequent repeats.

| VRS-Matcher workload | VRS-Matcher median (s) | VCFtools direct median (s) | Relative performance (% of VCFtools speed) |
|---|---:|---:|---:|
| Indexed all-pairs matching only | 0.157489 | 0.006001 | 3.8 |
| Index build + all-pairs matching | 0.320837 | 0.006001 | 1.9 |

The indexed all-pairs comparison excludes VCF reading for VRS-Matcher, while
the VCFtools measurement includes reading the VCF; it is not algorithm-only.
The second comparison includes VRS-Matcher index construction and is the
one-shot pipeline view. VCFtools reads the prepared VCF directly. Pair counts
also differ as noted above. The `query` stage is a separate one-versus-all
workload and has no single-query VCFtools counterpart in this benchmark.

Identity-only DB bytes: 61440.
Genotype + allele DB bytes: 61440.
Annotated VCF bytes: 1413.

See manifest.json, raw logs, pairs.tsv, and timings.tsv for reproducibility.
No biological accuracy claim is made by a smoke run.

The identical small database sizes reflect SQLite page allocation, not zero
storage cost for reference calls. The fixture is dominated by startup and is
far too small to estimate a cohort-scale crossover or memory growth. Downloads,
reference storage, preparation, and VRS annotation are excluded. VCFtools emits
ordered pairs and diagonals; the Python all-pairs worker computes each distinct
unordered pair once. Single-query timing is a separate workload.

## Raw timing observations

```tsv
repeat	stage	wall_seconds	cpu_seconds	peak_rss_bytes	exit_code
0	vcftools	0.01487416698364541	0.003306	2228224	0
0	load	0.16876425000373274	0.149387	39616512	0
0	query	0.16255191701930016	0.142382	37683200	0
0	all	0.15374716703081504	0.13974999999999999	37994496	0
1	load	0.1633483330369927	0.146571	39092224	0
1	query	0.1558726669754833	0.13869399999999998	37945344	0
1	all	0.15748908295063302	0.14199199999999998	37945344	0
1	vcftools	0.006063290988095105	0.003641	2080768	0
2	vcftools	0.004982999991625547	0.003181	2228224	0
2	load	0.15864533302374184	0.144585	39452672	0
2	query	0.1531172919785604	0.140009	38060032	0
2	all	0.15075987501768395	0.137435	37978112	0
3	load	0.16076916601741686	0.14354699999999998	39632896	0
3	query	0.15008325001690537	0.135546	37814272	0
3	all	0.20082254201406613	0.145462	37650432	0
3	vcftools	0.005788750015199184	0.0034609999999999997	2228224	0
4	vcftools	0.007960540999192744	0.00345	2228224	0
4	load	0.1662689580116421	0.14823999999999998	39157760	0
4	query	0.152376624988392	0.13817	37437440	0
4	all	0.15819916600594297	0.14196399999999998	37552128	0
5	load	0.16515049996087328	0.148378	39354368	0
5	query	0.14983995800139382	0.136247	37978112	0
5	all	0.1510751669993624	0.13701799999999997	37748736	0
5	vcftools	0.006001374975312501	0.003516	2080768	0
```

## Input fingerprints

- `tests/data/king/oracle.vcf`: `03e1ef69f2dcb7dba19dfccd4d9a1454fa25c71baa0594ccd4fbb7bbfc694b28`
- `tests/data/king/panel.tsv`: `3b3c93a47d9ceaab7e56ec4de52f79a987211510c53dfaf830dd79216b2e6aa4`
- `uv.lock`: `0e597f136c3dead769bc01239c03c73e9e84c428ead18333317fd73412b8ed17`

## Remaining biological work

Run the same tooling on a frozen genome-wide panel with verified independently
measured replicates, parent–offspring, siblings, and unrelated donors. Supply the
sample/donor/family/ancestry and relationship TSVs described in the
[usage guide](how-to-king-robust.md). The runner emits relationship and
ancestry/callable-overlap summaries, duplicate retrieval with ties, and
family-block bootstrap intervals. Those report paths are exercised with
synthetic labels in the test suite; synthetic labels do not satisfy biological
acceptance E. Record upstream annotation/preparation costs and evaluate larger
cohorts before drawing accuracy or end-to-end performance conclusions.
