# vrs-matcher

[![CI](https://github.com/ohsu-comp-bio/vrs-matcher/actions/workflows/ci.yml/badge.svg)](https://github.com/ohsu-comp-bio/vrs-matcher/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ohsu-comp-bio/vrs-matcher/graph/badge.svg)](https://codecov.io/gh/ohsu-comp-bio/vrs-matcher)

Match samples based on VRS identifiers.

## Requirements

- [uv](https://docs.astral.sh/uv/)
- Python 3.12+

## Quickstart

```bash
uv sync --all-groups
uv run pytest
uv run ruff check .
uv run ruff format .
```

## Project Layout

```text
src/vrs_matcher/    # package code
tests/              # test suite
pyproject.toml      # project metadata + tool configuration
```
