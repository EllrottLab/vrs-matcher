# vrs-matcher

[![CI](https://github.com/EllrottLab/vrs-matcher/actions/workflows/ci.yml/badge.svg)](https://github.com/EllrottLab/vrs-matcher/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/EllrottLab/vrs-matcher/graph/badge.svg)](https://codecov.io/gh/EllrottLab/vrs-matcher)

Match samples based on [GA4GH VRS](https://www.ga4gh.org/product/variation-representation/) identifiers.

## Quickstart

```bash
git clone https://github.com/EllrottLab/vrs-matcher

cd vrs-matcher

uv run vrs-matcher
```


## Examples

| VCF                       | Notebook                               |
|---------------------------|----------------------------------------|
| [example-cohort.vcf][vcf] | [![Open in Colab][colab-badge]][colab] |

[vcf]: ./examples/example-cohort.vcf
[colab-badge]: https://colab.research.google.com/assets/colab-badge.svg
[colab]: https://colab.research.google.com/github/EllrottLab/vrs-matcher/blob/8a620f6/examples/vrs-matcher.ipynb

- [Sample identity confirmation with your own VCFs](docs/how-to-sample-identity-confirmation.md)
- [Cohort deduplication and data release QC with your own VCFs](docs/how-to-cohort-dedup-qc.md)
- [KING-robust kinship estimation](docs/how-to-king-robust.md) — opt-in indexing of
  called reference and alternate SNP genotypes, pairwise estimates, and cohort ranking.

## Development

### Requirements

- [uv](https://docs.astral.sh/uv/)
- [Python 3.13+](https://www.python.org/downloads/)

### Tests

```bash
uv run pytest

uv run ruff check .

uv run ruff format .
```

## Plugin algorithms

`vrs-matcher` supports pluggable matching algorithms.

- Use a built-in plugin such as `identity` when the default workflow is enough.
- Use `king-robust` for kinship estimates after loading with
  `--index-genotypes --panel PANEL.tsv`; allele-only indexes require re-ingestion.
- Use a local script plugin when you want to prototype a lab- or study-specific matcher.
- Use an entry-point plugin when you want to distribute a reusable matcher as a Python package.

Plugins change the **matching and ranking logic**, not the ingestion pipeline.
Your VCFs are still loaded into the same SQLite-backed sample index; the plugin
controls how indexed samples are compared after loading.

```bash
uv run vrs-matcher plugins list
uv run vrs-matcher load-samples examples/example-cohort.vcf --db matches.db
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --algorithm identity
uv run vrs-matcher match-sample SAMPLE_A --db matches.db --plugin-file examples/plugins/jaccard_floor_plugin.py
```

See `docs/plugins.md` for:

- when to write a plugin,
- a minimal copy-pasteable plugin template,
- example bioinformatics plugin ideas such as rare-variant-weighted identity confirmation and candidate-gene-only sample matching,
- how to test a plugin on known samples, and
- how to package a plugin for reuse.

## Integration tests (opt-in)

Integration tests are marked `integration` and are skipped by default. Install
the integration dependencies before running them:

```bash
uv sync --group dev --group integration
```

Available integration tests:

- **1000 Genomes population signal** (`tests/integration/test_1kg_cluster.py`):
  downloads a chr22 VCF slice and annotates it with VRS IDs. It requires
  internet access and a local
  [seqrepo](https://github.com/biocommons/biocommons.seqrepo) data instance
  (one-time download, ~10 GB):

  ```bash
  # one-time seqrepo download (skip if already present)
  scripts/setup_integration_data.sh

  export GA4GH_VRS_DATAPROXY_URI=seqrepo+file://$HOME/.local/share/seqrepo/2024-12-20
  uv run pytest tests/integration/test_1kg_cluster.py --run-integration
  ```

- **KING/VCFtools comparison** (`tests/integration/test_king_vcftools.py`):
  compares KING-robust output with VCFtools `--relatedness2` using the small
  synthetic fixture in `tests/data/king/`. Install VCFtools 0.1.16 or newer
  using your system package manager:

  ```bash
  # macOS (Homebrew)
  brew install vcftools

  # Ubuntu/Debian
  sudo apt-get update
  sudo apt-get install vcftools
  ```

  Confirm the installed version with `vcftools --version`. If the executable
  is not on `PATH`, set `VCFTOOLS` to its path:

  ```bash
  VCFTOOLS=/opt/homebrew/bin/vcftools uv run pytest \
    tests/integration/test_king_vcftools.py --run-integration
  ```

- **Authorized GREGoR phenotype bundle** (`tests/integration/test_gregor.py`):
  validates phenotype import, observation enrichment, and cohort selection
  against already staged local data. This test does not download GREGoR data;
  only run it in an approved environment with authorized inputs. Set:

  ```bash
  export VRS_MATCHER_GREGOR_INTEGRATION=1
  export VRS_MATCHER_GREGOR_MANIFEST="/approved/path/manifest.json"
  export VRS_MATCHER_GREGOR_DB="/approved/path/genetic-index.db"
  export VRS_MATCHER_GREGOR_EXPECTED="/approved/path/expected-results.json"

  uv run pytest tests/integration/test_gregor.py --run-integration
  ```

  `VRS_MATCHER_GREGOR_MANIFEST` points to the bundle JSON manifest; file paths
  in it are resolved relative to the manifest's directory. `DB` is the verified
  genetic index referenced by the manifest, and `EXPECTED` is JSON containing
  the expected `sample_ids`, `observations`, `counts`, `predicate`, and
  `cohort_sample_ids`. Enabling the test without valid paths or data fails
  rather than silently skipping. Keep protected data and manifests in approved
  storage; do not commit them or credentials.

To run every integration test, configure seqrepo as above and make VCFtools
available, then run:

```bash
RUN_INTEGRATION_TESTS=1 uv run pytest -m integration --run-integration
```

To use a different seqrepo location/version, set `SEQREPO_ROOT` and
`SEQREPO_INSTANCE` before running `scripts/setup_integration_data.sh`, then
export the matching `GA4GH_VRS_DATAPROXY_URI`.

### Integration troubleshooting

- `Set GA4GH_VRS_DATAPROXY_URI to run integration tests`:
  export `GA4GH_VRS_DATAPROXY_URI` before pytest.
- `Could not initialise SeqRepo data proxy`:
  verify the directory in `GA4GH_VRS_DATAPROXY_URI` exists and contains the
  seqrepo snapshot you pulled.
- `Unable to fetch 1KGP remote VCF slice`:
  check internet access to `ftp.1000genomes.ebi.ac.uk` and retry.

The integration test downloads a chr22 region from 1000 Genomes, annotates it
with VRS IDs via `ga4gh.vrs`, ingests it into SQLite, and checks that mean
intra-super-population Jaccard is higher than inter-super-population Jaccard.

## Project Layout

```text
./vrs-matcher
├── pyproject.toml   # project metadata + tool configuration
├── src/             # package code
├── tests/           # test suite
├── docs/            # usage guides (e.g. plugins.md)
├── examples/        # runnable examples + sample data (see examples/README.md)
└── scripts/         # helper / maintenance scripts
```
