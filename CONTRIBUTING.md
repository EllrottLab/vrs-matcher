# Contributing

Thank you for your interest in contributing to vrs-matcher!

## Development setup

Install [uv](https://docs.astral.sh/uv/), then:

```bash
git clone https://github.com/ohsu-comp-bio/vrs-matcher.git
cd vrs-matcher
uv sync --all-groups
```

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
```

## Pull requests

- Keep PRs focused on a single change
- Include tests for new behaviour
- Ensure all checks pass
