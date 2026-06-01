"""Notebook smoke tests for the bundled example notebook(s)."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

from vrs_matcher.db import list_samples, open_db
from vrs_matcher.loader import load_samples
from vrs_matcher.matcher import match_against_all, match_pair


def _notebook_code_sources(path: Path) -> list[str]:
    """Return the source of all code cells in a notebook file."""

    notebook = json.loads(path.read_text(encoding="utf-8"))
    return ["".join(cell.get("source", [])) for cell in notebook.get("cells", []) if cell.get("cell_type") == "code"]


def test_example_notebook_contains_matcher_workflow() -> None:
    """Verify the example notebook still documents the expected matcher API."""

    notebook_path = Path("examples/vrs-matcher.ipynb")
    assert notebook_path.exists()

    sources = _notebook_code_sources(notebook_path)
    notebook_text = "\n".join(sources)

    assert "from vrs_matcher.loader import load_samples" in notebook_text
    assert "from vrs_matcher.matcher import match_pair" in notebook_text
    assert "from vrs_matcher.matcher import match_against_all" in notebook_text


def test_example_notebook_matching_smoke(tmp_path: Path) -> None:
    """Execute the notebook's core matching workflow against the bundled example cohort."""

    # The notebook downloads the same example cohort files into its working
    # directory; this test uses the checked-in copies so it can run offline.
    vcf_path = Path("examples/example-cohort.vcf.gz")
    assert vcf_path.exists()

    with TemporaryDirectory(dir=tmp_path) as workdir:
        db_path = Path(workdir) / "cohort.db"
        load_samples(str(vcf_path), str(db_path))

        conn = open_db(db_path)
        try:
            samples = list_samples(conn)
            assert len(samples) >= 2

            pair = match_pair(conn, "SAMPLE_A", "SAMPLE_B")
            assert pair.jaccard == 2 / 3
            assert pair.weighted_concordance == 1.0

            ranked = match_against_all(conn, "SAMPLE_A")
            assert ranked
            assert ranked[0].sample_b == "SAMPLE_B"
        finally:
            conn.close()

