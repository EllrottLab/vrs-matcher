# Contributing

Thank you for your interest in contributing to vrs-matcher!

## Development setup

Install [uv](https://docs.astral.sh/uv/), then:

```bash
git clone https://github.com/ohsu-comp-bio/vrs-matcher.git
cd vrs-matcher
uv sync --all-groups
```

If you plan to run the real-data integration test, install the integration
group and a local [seqrepo](https://github.com/biocommons/biocommons.seqrepo)
data instance:

```bash
uv sync --group dev --group integration

# one-time seqrepo download (~10 GB, cached across runs)
scripts/setup_integration_data.sh

# required by tests/integration/conftest.py
export GA4GH_VRS_DATAPROXY_URI=seqrepo+file://$HOME/.local/share/seqrepo/2024-12-20
```

To use a different snapshot or root directory, set `SEQREPO_INSTANCE` and/or
`SEQREPO_ROOT` before running `scripts/setup_integration_data.sh`.

## Making changes

1. Create a branch: `git checkout -b your-feature`
2. Make your changes in `src/vrs_matcher/`
3. Add or update tests in `tests/`
4. Run the checks below and fix any failures before opening a PR

## Checks

```bash
# Lint and format
uv run ruff check .
uv run ruff format .

# Tests
uv run pytest

# Optional integration test (network + local seqrepo required)
GA4GH_VRS_DATAPROXY_URI=seqrepo+file://$HOME/.local/share/seqrepo/2024-12-20 \
  RUN_INTEGRATION_TESTS=1 uv run pytest -m integration --run-integration
```

## Pull requests

- Keep PRs focused on a single change
- Include tests for new behaviour
- Ensure all checks pass
