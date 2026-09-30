"""KING acceptance A–C: real VCF parsing, atomic ingestion, and independent oracles."""

import json
from dataclasses import asdict
from pathlib import Path

import pytest
from click.testing import CliRunner

from vrs_matcher.cli import cli
from vrs_matcher.db import get_vrs_ids, insert_alleles, open_db, register_samples
from vrs_matcher.genotypes import read_panel
from vrs_matcher.loader import load_samples
from vrs_matcher.matcher import match_against_all, match_pair
from vrs_matcher.plugins import PluginContext, PluginError, list_plugins

DATA = Path(__file__).parent / "data" / "king"


def inputs(tmp_path, pairs=None):
    vcf = tmp_path / "calls.vcf"
    panel = tmp_path / "panel.tsv"
    text = (DATA / "oracle.vcf").read_text()
    if pairs is not None:
        header = [s for s in text.splitlines() if s.startswith("#")]
        records = [s.split("\t") for s in text.splitlines() if not s.startswith("#")]
        rows = []
        for record, pair in zip(records, pairs, strict=False):
            record[9:] = [f"{gt}:20:10" for gt in pair]
            rows.append("\t".join(record))
        text = "\n".join(header + rows) + "\n"
    vcf.write_text(text)
    panel.write_text((DATA / "panel.tsv").read_text())
    return vcf, panel, tmp_path / "index.db"


def load(vcf, panel, db, **kwargs):
    return load_samples(vcf, db, index_genotypes=True, panel=panel, **kwargs)


def result(db, a="A", b="B", **kwargs):
    conn = open_db(db)
    try:
        return match_pair(conn, a, b, algorithm="king-robust", **kwargs)
    finally:
        conn.close()


def test_oracle_and_reopen(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    summary = load(vcf, panel, db)
    r = result(db)
    assert (r.n_common, r.het_a, r.het_b, r.het_both, r.opposite_hom) == (8, 4, 6, 3, 1)
    assert r.kinship == pytest.approx(0, abs=1e-12)
    assert r.kinship_within_family == pytest.approx(0.1, abs=1e-12)
    assert r.ibs0_fraction == 0.125
    assert r == result(db)
    assert summary["dosage_0"] == 5
    assert summary["dosage_1"] == 10
    assert summary["dosage_2"] == 1
    conn = open_db(db)
    assert len(PluginContext(conn).get_called_genotypes("A")) == 8
    assert len(get_vrs_ids(conn, "A")) == 4  # Reference rows never enter identity sets.
    observation = conn.execute("SELECT * FROM genotype_observation WHERE sample_id='A'").fetchone()
    assert json.loads(observation["summary"])["calls"]["dosage_0"] == 4
    assert len(json.loads(observation["provenance"])["input_sha256"]) == 64
    conn.close()
    assert "king-robust" in list_plugins()


@pytest.mark.parametrize(
    ("pairs", "m", "score", "within", "opposite"),
    [
        ([("0/1", "0/1"), ("1/1", "1/1")], 2, 0.5, 0.5, 0),
        ([("0/1", "0/1"), ("0/0", "1/1")], 2, -0.5, -0.5, 1),
        ([("./.", "./.")], 0, None, None, 0),
        ([("0/0", "0/0")], 1, None, None, 0),
        ([("0/0", "0/1")], 1, None, 0, 0),
        ([("0/1", "0/1"), ("0/1", "./.")], 1, 0.5, 0.5, 0),
        ([("1/1", "./1"), ("1/1", "./.")], 0, None, None, 0),
        ([("1|0", "0|1")], 1, 0.5, 0.5, 0),
        ([("1", "0/1"), ("0/1/1", "0/1"), ("0/2", "0/1")], 0, None, None, 0),
    ],
)
def test_statistical_and_callability_cases(tmp_path, pairs, m, score, within, opposite):
    vcf, panel, db = inputs(tmp_path, pairs)
    load(vcf, panel, db)
    r = result(db)
    assert (r.n_common, r.kinship, r.kinship_within_family, r.opposite_hom) == (
        m,
        score,
        within,
        opposite,
    )
    assert r.status == ("unscorable" if score is None else "ok")
    assert (r.reason is not None) == (score is None)
    assert (r.within_family_reason is not None) == (within is None)
    json.dumps(asdict(r), allow_nan=False)


def test_phase_sample_order_record_order_and_ref_alt_invariance(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    load(vcf, panel, db)
    expected = result(db)
    swapped = result(db, "B", "A")
    assert (swapped.het_a, swapped.het_b) == (6, 4)
    assert swapped.kinship == expected.kinship
    assert result(db, "A", "A").kinship == 0.5
    header, records = [], []
    for line in vcf.read_text().splitlines():
        if line.startswith("#"):
            header.append(line)
            continue
        cells = line.split("\t")
        cells[3:5] = ["C", "A"]
        for i in (9, 10):
            gt, rest = cells[i].split(":", 1)
            gt = gt.translate(str.maketrans("01/", "10|"))
            cells[i] = f"{gt}:{rest}"
        records.append("\t".join(cells))
    vcf.write_text("\n".join(header + records[::-1]) + "\n")
    panel.write_text(panel.read_text().replace("\tA\tC\t", "\tC\tA\t"))
    other = tmp_path / "reversed.db"
    load(vcf, panel, other)
    r = result(other)
    assert (r.kinship, r.kinship_within_family, r.opposite_hom) == (
        expected.kinship,
        expected.kinship_within_family,
        expected.opposite_hom,
    )


def test_restrictions_and_unknown_ids(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    load(vcf, panel, db)
    candidates = frozenset(f"ga4gh:VA.{i:032d}" for i in range(1, 4))
    r = result(db, candidate_vrs_ids=candidates)
    assert r.n_common == 3 and r.kinship == 0.5
    assert result(db, candidate_vrs_ids=frozenset()).n_common == 0
    assert result(db, candidate_vrs_ids=frozenset([f"ga4gh:VA.{8:032d}"])).opposite_hom == 1
    with pytest.raises(PluginError, match="belong"):
        result(db, candidate_vrs_ids=frozenset(["unknown"]))


def test_qc_and_unsupported_records(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    lines = vcf.read_text().splitlines()
    records = [s.split("\t") for s in lines if not s.startswith("#")]
    # Same QC policy for reference and alternate calls. Missing QC is allowed.
    records[0][9:] = ["0/0:19:10", "1/1:19:10"]
    records[1][9:] = ["0/0:20:9", "1/1:20:9"]
    records[2][9:] = ["0/0:.:.", "1/1:.:."]
    records[3][6] = "q10"
    records[4][4] = "CC"  # Unsupported indel.
    records[5][4] = "C,G"
    records[6][0] = "chrX"
    records[7][6] = "."
    vcf.write_text(
        "\n".join([s for s in lines if s.startswith("#")] + ["\t".join(s) for s in records]) + "\n"
    )
    summary = load(vcf, panel, db, dp_threshold=10)
    assert result(db).n_common == 2
    assert result(db).opposite_hom == 2
    assert summary["excluded_low_gq"] == summary["excluded_low_dp"] == 2
    assert summary["missing_gq"] == summary["missing_dp"] == 2
    assert summary["excluded_filter_records"] == 1
    assert summary["excluded_unsupported_records"] == 2
    assert summary["excluded_nonpanel_contig_records"] == 1
    assert summary["eligible_calls"] == 8
    assert summary["records"] == 8


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        ("VRS_Allele_IDs=ga4gh:VA." + "0" * 31 + "8", ".", "annotation"),
        ("VRS_Allele_IDs=ga4gh:VA." + "0" * 31 + "8", "VRS_Allele_IDs=.", "annotation"),
        ("VRS_Allele_IDs=ga4gh:VA." + "0" * 31 + "8", "VRS_Allele_IDs=x,y", "annotation"),
        ("\t8\t.\tA\tC\t", "\t8\t.\tC\tA\t", "REF/ALT"),
        ("##reference=synthetic-test-reference", "##reference=other", "reference"),
        ("refget=ga4gh:SQ." + "A" * 32, "refget=ga4gh:SQ." + "B" * 32, "refget"),
    ],
)
def test_failed_load_rolls_back(tmp_path, old, new, message):
    vcf, panel, db = inputs(tmp_path)
    vcf.write_text(vcf.read_text().replace(old, new))
    with pytest.raises(ValueError, match=message):
        load(vcf, panel, db)
    conn = open_db(db)
    for table in (
        "samples",
        "sample_allele",
        "genotype_call",
        "genotype_observation",
        "genotype_panel",
    ):
        assert conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0] == 0
    conn.close()


def test_duplicate_records_and_append_conflicts(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    original = vcf.read_text()
    vcf.write_text(original + original.splitlines()[-1] + "\n")
    with pytest.raises(ValueError, match="Duplicate"):
        load(vcf, panel, db)
    vcf.write_text(original)
    load(vcf, panel, db)
    before = result(db)
    with pytest.raises(ValueError, match="already exists"):
        load(vcf, panel, db)
    vcf.write_text(original.replace("\tA\tB\n", "\tC\tD\n"))
    with pytest.raises(ValueError, match="QC policy"):
        load(vcf, panel, db, gq_threshold=30)
    assert before == result(db)
    load(vcf, panel, db)
    assert result(db, "A", "C").kinship == 0.5


def test_injected_failure_during_identity_projection(tmp_path, monkeypatch):
    from vrs_matcher import genotypes

    vcf, panel, db = inputs(tmp_path)

    def fail(*args, **kwargs):
        raise RuntimeError("injected failure")

    monkeypatch.setattr(genotypes, "insert_alleles", fail)
    with pytest.raises(RuntimeError, match="injected"):
        load(vcf, panel, db)
    conn = open_db(db)
    assert conn.execute("SELECT COUNT(*) FROM genotype_call").fetchone()[0] == 0
    assert conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 0
    conn.close()


def test_legacy_and_mixed_indexes(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    load_samples(vcf, db)
    with pytest.raises(PluginError, match="re-ingestion"):
        result(db)
    with pytest.raises(ValueError, match="already exists"):
        load(vcf, panel, db)
    fresh = tmp_path / "fresh.db"
    load(vcf, panel, fresh)
    conn = open_db(fresh)
    register_samples(conn, ["LEGACY"])
    assert match_pair(conn, "A", "LEGACY").jaccard == 0
    with pytest.raises(PluginError, match="re-ingestion"):
        match_against_all(conn, "A", algorithm="king-robust")
    with pytest.raises(ValueError, match="overwrite"):
        insert_alleles(conn, [("A", "x", "0/1", "HET", "chr1", 1, None, None, None)])
    conn.close()


def test_ranking_cli_and_json(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    runner = CliRunner()
    loaded = runner.invoke(
        cli, ["load-samples", str(vcf), "--db", str(db), "--index-genotypes", "--panel", str(panel)]
    )
    assert loaded.exit_code == 0, loaded.output
    original = vcf.read_text()
    vcf.write_text(original.replace("\tA\tB\n", "\tC\tD\n"))
    load(vcf, panel, db)
    vcf.write_text(original.replace("\tA\tB\n", "\tEMPTY\tREF\n"))
    lines = vcf.read_text().splitlines()
    for i, line in enumerate(lines):
        if not line.startswith("#"):
            cells = line.split("\t")
            cells[9:] = ["./.:20:10", "0/0:20:10"]
            lines[i] = "\t".join(cells)
    vcf.write_text("\n".join(lines) + "\n")
    load(vcf, panel, db)
    conn = open_db(db)
    ranking = match_against_all(conn, "A", algorithm="king-robust", top_n=2)
    assert [r.sample_b for r in ranking.matches] == ["C", "B"]
    assert [r.sample_b for r in ranking.unscorable] == ["EMPTY", "REF"]
    for r in ranking.matches + ranking.unscorable:
        assert r == match_pair(conn, "A", r.sample_b, algorithm="king-robust")
    conn.close()
    pair_args = ["match-samples", "A", "B", "--db", str(db), "--algorithm", "king-robust"]
    pair = runner.invoke(cli, pair_args + ["--json"])
    assert pair.exit_code == 0, pair.output
    assert json.loads(pair.output)["kinship_within_family"] == 0.1
    assert "jaccard" not in pair.output
    assert "kinship:" in runner.invoke(cli, pair_args).output
    search = runner.invoke(
        cli,
        [
            "match-sample",
            "A",
            "--db",
            str(db),
            "--algorithm",
            "king-robust",
            "--top",
            "1",
            "--json",
        ],
    )
    assert search.exit_code == 0, search.output
    assert len(json.loads(search.output)["matches"]) == 1
    assert len(json.loads(search.output)["unscorable"]) == 2
    rejected = runner.invoke(
        cli, ["shared-variants", "A", "B", "--db", str(db), "--algorithm", "king-robust"]
    )
    assert rejected.exit_code != 0 and "requires an allele" in rejected.output
    legacy = runner.invoke(cli, ["match-samples", "A", "B", "--db", str(db), "--json"])
    assert legacy.exit_code == 0 and "jaccard" in json.loads(legacy.output)


def test_batched_and_streamed_king_pairs_and_top_n_ties(tmp_path):
    from itertools import combinations

    from vrs_matcher.king import KingRobustPlugin
    from vrs_matcher.plugins import PluginContext

    vcf, panel, db = inputs(tmp_path)
    lines = vcf.read_text().splitlines()
    changed = []
    for line in lines:
        if line.startswith("#CHROM"):
            changed.append("\t".join(line.split("\t")[:9] + ["A", "B", "C", "D"]))
        elif not line.startswith("#"):
            cells = line.split("\t")
            cells[9:] = [cells[9]] * 4
            changed.append("\t".join(cells))
        else:
            changed.append(line)
    vcf.write_text("\n".join(changed) + "\n")
    load(vcf, panel, db)

    conn = open_db(db)
    context = PluginContext(conn)
    batched = context.get_called_genotypes_many(["A", "B", "A"], batch_size=1)
    assert list(batched) == ["A", "B"]
    assert batched["A"] == context.get_called_genotypes("A")

    plugin = KingRobustPlugin()
    all_pairs = list(plugin.iter_all_pairs(context, batch_size=2))
    assert [(r.sample_a, r.sample_b) for r in all_pairs] == list(combinations("ABCD", 2))
    assert all(
        pair == match_pair(conn, pair.sample_a, pair.sample_b, algorithm="king-robust")
        for pair in all_pairs
    )
    complete = match_against_all(conn, "A", algorithm="king-robust")
    top_two = match_against_all(conn, "A", algorithm="king-robust", top_n=2)
    assert (
        [r.sample_b for r in top_two.matches]
        == [r.sample_b for r in complete.matches[:2]]
        == ["B", "C"]
    )
    conn.close()


@pytest.mark.parametrize("mutation", ["duplicate", "empty", "sequence", "reference", "nonauto"])
def test_panel_validation(tmp_path, mutation):
    _, panel, _ = inputs(tmp_path)
    text = panel.read_text()
    if mutation == "duplicate":
        text += text.splitlines()[-1] + "\n"
    elif mutation == "empty":
        text = text.splitlines()[0] + "\n"
    elif mutation == "sequence":
        text = text.replace("ga4gh:SQ." + "A" * 32, "chr1", 1)
    elif mutation == "reference":
        text = text.replace("synthetic-test-reference", "different", 1)
    else:
        text = text.replace("\t1\tga4gh:", "\tX\tga4gh:", 1)
    panel.write_text(text)
    with pytest.raises(ValueError):
        read_panel(panel)


def test_absent_marker_and_aliases(tmp_path):
    vcf, panel, db = inputs(tmp_path, [("0/1", "0/1"), ("0/0", "1/1")])
    original = vcf.read_text()
    # Real alias resolution requires matching sequence digest, not just 'chr' stripping.
    vcf.write_text(original.replace("chr1", "1"))
    load(vcf, panel, db)
    vcf.write_text("\n".join(original.splitlines()[:-1]).replace("\tA\tB\n", "\tC\tD\n") + "\n")
    load(vcf, panel, db)
    r = result(db, "B", "C")
    assert r.n_common == 1 and r.opposite_hom == 0


def test_required_flags_and_threshold_validation(tmp_path):
    vcf, panel, db = inputs(tmp_path)
    with pytest.raises(ValueError, match="requires --panel"):
        load_samples(vcf, db, index_genotypes=True)
    with pytest.raises(ValueError, match="requires --index-genotypes"):
        load_samples(vcf, db, panel=panel)
    with pytest.raises(ValueError, match="finite"):
        load(vcf, panel, db, gq_threshold=float("nan"))
    with pytest.raises(ValueError, match="legacy filters"):
        load(vcf, panel, db, candidate_vrs_ids=set())


def test_benchmark_preparation_and_label_report(tmp_path):
    """Exercise the report's external-label accounting without pretending labels are biological."""
    import importlib.util
    from argparse import Namespace

    spec = importlib.util.spec_from_file_location(
        "benchmark_king", Path("scripts/benchmark_king.py")
    )
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    vcf, panel, _ = inputs(tmp_path)
    args = Namespace(
        vcf=vcf,
        panel=panel,
        samples=tmp_path / "samples.tsv",
        relationships=tmp_path / "relationships.tsv",
        output=tmp_path,
    )
    benchmark.validate_prepared_input(args)
    vcf.write_text(vcf.read_text().replace(":20:10", ":19:10"))
    with pytest.raises(ValueError, match="QC masking"):
        benchmark.validate_prepared_input(args)
    args.samples.write_text(
        "sample\tdonor\tfamily\tancestry\n"
        "A\td1\tf1\ttest\nB\td1\tf1\ttest\n"
        "C\td2\tf1\ttest\nD\td3\tf1\ttest\n"
        "E\td4\tf2\ttest\nF\td4\tf2\ttest\n"
    )
    args.relationships.write_text(
        "sample_a\tsample_b\trelationship\n"
        "A\tB\tduplicate\nA\tC\tparent_offspring\n"
        "C\tD\tsiblings\nA\tE\tunrelated\nE\tF\tduplicate\n"
    )
    import itertools

    scores = [
        {
            "sample_a": a,
            "sample_b": b,
            "n_common": "8",
            "kinship": ".5" if a + b in ("AB", "EF") else ".1",
            "kinship_within_family": ".5" if a + b in ("AB", "EF") else ".1",
            "vcftools_phi": ".5" if a + b in ("AB", "EF") else ".1",
        }
        for a, b in itertools.combinations("ABCDEF", 2)
    ]
    benchmark.summarize_labels(args, scores)
    assert len(benchmark.read_tsv(tmp_path / "labeled-pairs.tsv")) == 5
    retrieved = benchmark.read_tsv(tmp_path / "duplicate-retrieval.tsv")
    assert {row["metric"] for row in retrieved} == {
        "kinship",
        "kinship_within_family",
        "vcftools_phi",
    }
    assert all(r["unique_correct"] == "1" for r in retrieved)
    assert json.loads((tmp_path / "uncertainty.json").read_text())["kinship"][
        "accuracy_95_percent_interval"
    ] == [1, 1]


def test_benchmark_workload_summaries_use_matched_repeats():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "benchmark_king", Path("scripts/benchmark_king.py")
    )
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    timings = [
        {"repeat": 0, "stage": "load", "wall_seconds": 20},
        {"repeat": 0, "stage": "all", "wall_seconds": 30},
        {"repeat": 0, "stage": "vcftools", "wall_seconds": 2},
        {"repeat": 1, "stage": "load", "wall_seconds": 2},
        {"repeat": 1, "stage": "all", "wall_seconds": 3},
        {"repeat": 1, "stage": "vcftools", "wall_seconds": 1},
        {"repeat": 2, "stage": "load", "wall_seconds": 8},
        {"repeat": 2, "stage": "all", "wall_seconds": 2},
        {"repeat": 2, "stage": "vcftools", "wall_seconds": 1},
    ]
    summary = benchmark.summarize_workloads(timings)
    pipeline = summary["index_plus_all_pairs"]
    indexed = summary["indexed_all_pairs"]
    assert pipeline["plugin_median_seconds"] == 7.5
    assert pipeline["vcftools_median_seconds"] == 1
    assert pipeline["relative_performance_percent"] == pytest.approx(100 / 7.5)
    assert pipeline["basis"] == "subsequent runs"
    assert indexed["plugin_median_seconds"] == 2.5
    assert indexed["relative_performance_percent"] == pytest.approx(40)
    assert summary["index_build_median_seconds"] == 5
    assert summary["query_median_seconds"] is None
    first_only = [row for row in timings if row["repeat"] == 0]
    first_summary = benchmark.summarize_workloads(first_only)
    assert first_summary["index_plus_all_pairs"]["basis"] == "first run only"
    assert first_summary["indexed_all_pairs"]["plugin_median_seconds"] == 30
