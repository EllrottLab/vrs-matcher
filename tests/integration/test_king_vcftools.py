"""Acceptance D: opt-in comparison to a pinned external VCFtools executable."""

import csv
import gzip
import os
import re
import shutil
import subprocess
from decimal import Decimal
from pathlib import Path

import pytest

from vrs_matcher.db import open_db
from vrs_matcher.loader import load_samples
from vrs_matcher.matcher import match_pair

DATA = Path(__file__).parents[1] / "data" / "king"


@pytest.mark.integration
@pytest.mark.parametrize("missing", [False, True])
def test_vcftools_relatedness2(tmp_path, missing):
    executable = os.environ.get("VCFTOOLS") or shutil.which("vcftools")
    if not executable:
        pytest.skip("Install VCFtools 0.1.16 or newer, or set VCFTOOLS to its executable")
    version = subprocess.run([executable, "--version"], capture_output=True, text=True, check=True)
    version_match = re.fullmatch(r"VCFtools \((\d+)\.(\d+)\.(\d+)\)", version.stdout.strip())
    assert version_match is not None, f"Unrecognized VCFtools version: {version.stdout.strip()!r}"
    assert tuple(map(int, version_match.groups())) >= (0, 1, 16), (
        f"VCFtools 0.1.16 or newer is required; found {version.stdout.strip()}"
    )
    text = (DATA / "oracle.vcf").read_text()
    if missing:
        header = [s for s in text.splitlines() if s.startswith("#")]
        rows = [s.split("\t") for s in text.splitlines() if not s.startswith("#")][:2]
        rows[1][10] = "./.:20:10"
        text = "\n".join(header + ["\t".join(r) for r in rows]) + "\n"
    vcf = tmp_path / "panel.vcf.gz"
    with gzip.open(vcf, "wt") as out:
        out.write(text)
    db = tmp_path / "king.db"
    load_samples(vcf, db, index_genotypes=True, panel=DATA / "panel.tsv")
    baseline = subprocess.run(
        [executable, "--gzvcf", str(vcf), "--relatedness2", "--out", str(tmp_path / "baseline")],
        capture_output=True,
        text=True,
        check=True,
    )
    (tmp_path / "command.log").write_text(version.stdout + baseline.stdout + baseline.stderr)
    conn = open_db(db)
    try:
        with (tmp_path / "baseline.relatedness2").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        assert len(rows) == 4  # Includes self comparisons and both pair directions.
        for row in rows:
            r = match_pair(conn, row["INDV1"], row["INDV2"], algorithm="king-robust")
            assert int(row["N_AaAa"]) == r.het_both
            assert int(row["N_AAaa"]) == r.opposite_hom
            raw = row["RELATEDNESS_PHI"]
            tolerance = 0.5 * float(Decimal(10) ** Decimal(raw).as_tuple().exponent) + 1e-12
            if missing and row["INDV1"] != row["INDV2"]:
                assert r.kinship == r.kinship_within_family == 0.5
                assert float(raw) == pytest.approx(1 / 3, abs=tolerance)
                assert sorted([int(row["N1_Aa"]), int(row["N2_Aa"])]) == [1, 2]
            else:
                assert int(row["N1_Aa"]) == r.het_a
                assert int(row["N2_Aa"]) == r.het_b
                assert float(raw) == pytest.approx(r.kinship_within_family, abs=tolerance)
    finally:
        conn.close()
