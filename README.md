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
[colab]: https://colab.research.google.com/github/EllrottLab/vrs-matcher/blob/be7e77e/examples/vrs-matcher.ipynb

## Development

### Requirements

- [uv](https://docs.astral.sh/uv/)
- [Python 3.12+](https://www.python.org/downloads/)

### Tests

```bash
uv run pytest

uv run ruff check .

uv run ruff format .
```

## Integration test (opt-in)

The 1000 Genomes end-to-end test is marked `integration` and is skipped by
default. It requires a local [seqrepo](https://github.com/biocommons/biocommons.seqrepo)
data instance (one-time download, ~10 GB) and internet access for the 1KGP VCF.

```bash
# install only what the integration test needs
uv sync --group dev --group integration

# one-time seqrepo download (skip if already present)
scripts/setup_integration_data.sh

# run just integration tests
export GA4GH_VRS_DATAPROXY_URI=seqrepo+file://$HOME/.local/share/seqrepo/2024-12-20
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
├── src
│   └── vrs_matcher  # package code
└── tests            # test suite
```
