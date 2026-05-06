# vrs-matcher

[![CI](https://github.com/EllrottLab/vrs-matcher/actions/workflows/ci.yml/badge.svg)](https://github.com/EllrottLab/vrs-matcher/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/EllrottLab/vrs-matcher/graph/badge.svg)](https://codecov.io/gh/EllrottLab/vrs-matcher)

Match samples based on [GA4GH VRS](https://www.ga4gh.org/product/variation-representation/) identifiers.

## Quickstart

```bash
git clone https://github.com/EllrottLab/vrs-matcher

cd vrs-matcher

uv sync

uv run vrs-matcher
```

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

### Project Layout

```text
./vrs-matcher
├── pyproject.toml   # project metadata + tool configuration
├── src
│   └── vrs_matcher  # package code
└── tests            # test suite
```
