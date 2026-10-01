"""Synthetic clinical contracts and genetic parity; no network or clinical data."""

import json
import shutil
import sqlite3
from contextlib import closing
from dataclasses import asdict
from pathlib import Path

import pytest
from click.testing import CliRunner

from vrs_matcher.cli import cli
from vrs_matcher.genotypes import file_sha256, load_genotypes
from vrs_matcher.matcher import match_against_all, match_pair
from vrs_matcher.phenotypes import (
    Phenotypes,
    export_pairs,
    import_snapshot,
    match_report,
    read_only,
    snapshot_diff,
)

FIXTURE = Path(__file__).parent / "data/gregor"


@pytest.fixture
def bundle(tmp_path):
    root = tmp_path / "inputs"
    shutil.copytree(FIXTURE, root)
    db = tmp_path / "genetic.db"
    load_genotypes(root / "calls.vcf", db, root / "panel.tsv")
    manifest = root / "bundle.json"
    data = json.loads(manifest.read_text())
    data["genetic"]["sha256"] = file_sha256(db)
    manifest.write_text(json.dumps(data))
    return db, tmp_path / "clinical.db", manifest


def refresh(manifest):
    data = json.loads(manifest.read_text())
    for source in data["sources"]:
        for key in ("participant", "phenotype", "mapping", "reconciliation"):
            if key in source:
                source[key]["sha256"] = file_sha256(manifest.parent / source[key]["path"])
    manifest.write_text(json.dumps(data))


def append(path, text):
    with path.open("a") as stream:
        stream.write(text)


@pytest.mark.parametrize("version", ["1.11", "1.12"])
def test_import_enrichment_presence_and_replay(bundle, version):
    db, sidecar, manifest = bundle
    data = json.loads(manifest.read_text())
    data["model_version"] = version
    manifest.write_text(json.dumps(data))
    before = file_sha256(db)
    result = import_snapshot(*bundle)
    assert result["phenotype_participant"] == 6
    assert result["unmapped_observations"] == 1
    assert result["presence"] == {"Absent": 2, "Present": 2, "Unknown": 1}
    assert import_snapshot(*bundle)["already_present"]
    assert file_sha256(db) == before
    with closing(Phenotypes(sidecar, result["snapshot_id"], db)) as reader:
        report = reader.enrich(["S0", "S1", "S2", "S5"])
        assert report["counts"]["people"] == 2
        assert (
            report["observations"]["S0"]["person_key"] == report["observations"]["S1"]["person_key"]
        )
        assert report["observations"]["S5"]["status"] == "unmapped"
        assert reader.cohort("HPO", "HP:0000001")["sample_ids"] == ["S0", "S1", "S3"]
        assert reader.cohort("HPO", "HP:0000001")["groups"]["unrecorded"] == 1
        assert reader.cohort("HPO", "HP:0000002")["sample_ids"] == []
        assert reader.cohort("MONDO", "HP:0000001")["sample_ids"] == []


def test_conflict_reconciliation_terra_and_snapshot_diff(bundle):
    db, sidecar, manifest = bundle
    first = import_snapshot(*bundle)["snapshot_id"]
    root = manifest.parent
    # Another source with the same local ID is independent until explicitly reconciled.
    (root / "other.tsv").write_text("entity:participant_id\tfamily_id\nP0\tF0\n")
    (root / "other_phen.tsv").write_text(
        "entity:phenotype_id\tparticipant_id\tterm_id\tpresence\tontology\n"
        'export-id\t{"entityType":"participant","entityName":"P0"}\tHP:0000001\tAbsent\tHPO\n'
    )
    (root / "empty.tsv").write_text("sample_id\tparticipant_id\n")
    data = json.loads(manifest.read_text())
    source = {
        "source_id": "second",
        "participant": {"path": "other.tsv"},
        "phenotype": {"path": "other_phen.tsv"},
        "mapping": {"path": "empty.tsv"},
    }
    data["sources"].append(source)
    manifest.write_text(json.dumps(data))
    refresh(manifest)
    independent = import_snapshot(*bundle)["snapshot_id"]
    with closing(Phenotypes(sidecar, independent)) as reader:
        assert not reader.enrich(["S0"])["participants"].popitem()[1]["terms"][0]["conflict"]
    (root / "reconcile.tsv").write_text("participant_id\tperson_key\nP0\tLINK\n")
    for s in data["sources"]:
        s["reconciliation"] = {"path": "reconcile.tsv"}
    manifest.write_text(json.dumps(data))
    refresh(manifest)
    second = import_snapshot(*bundle)["snapshot_id"]
    with closing(Phenotypes(sidecar, second)) as reader:
        p = reader.enrich(["S0"])["participants"].popitem()[1]
        assert len(p["records"]) == 2 and len(p["assertions"]) == 2
        assert p["terms"][0]["statuses"] == ["Absent", "Present"]
        assert reader.cohort("HPO", "HP:0000001", exclude_conflicts=True)["sample_ids"] == ["S3"]
    assert snapshot_diff(sidecar, first, second)["changes"]["assertions"]["added"]
    lines = (root / "phenotype.tsv").read_text().splitlines()
    (root / "phenotype.tsv").write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n")
    refresh(manifest)
    third = import_snapshot(*bundle)["snapshot_id"]
    assert third != second
    assert all(
        not v
        for change in snapshot_diff(sidecar, second, third)["changes"].values()
        for v in change.values()
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "presence",
        "namespace",
        "duplicate_mapping",
        "unknown_mapping",
        "reference",
        "schema",
        "checksum",
        "missing",
        "duplicate_header",
    ],
)
def test_bad_import_rolls_back_preserving_old_snapshot(bundle, mutation):
    _, sidecar, manifest = bundle
    old = import_snapshot(*bundle)["snapshot_id"]
    root = manifest.parent
    if mutation == "presence":
        append(root / "phenotype.tsv", "P0\tHP:0000001\tMaybe\tHPO\tx\n")
    elif mutation == "namespace":
        append(root / "phenotype.tsv", "P0\tMONDO:0000001\tPresent\tHPO\tx\n")
    elif mutation == "reference":
        append(root / "phenotype.tsv", "UNKNOWN\tHP:0000001\tPresent\tHPO\tx\n")
    elif mutation in ("duplicate_mapping", "unknown_mapping"):
        append(
            root / "mapping.tsv",
            ("S0" if mutation == "duplicate_mapping" else "UNKNOWN") + "\tP1\n",
        )
    elif mutation == "schema":
        data = json.loads(manifest.read_text())
        data["model_version"] = "9"
        manifest.write_text(json.dumps(data))
    elif mutation == "missing":
        append(root / "phenotype.tsv", "P0\t\tPresent\tHPO\tx\n")
    elif mutation == "duplicate_header":
        (root / "mapping.tsv").write_text("sample_id\tparticipant_id\tsample_id\n")
    else:
        append(root / "phenotype.tsv", "P0\tHP:0000001\tAbsent\tHPO\tx\n")
    if mutation != "checksum":
        refresh(manifest)
    with pytest.raises((ValueError, sqlite3.IntegrityError)):
        import_snapshot(*bundle)
    with closing(read_only(sidecar)) as conn:
        assert [r[0] for r in conn.execute("SELECT snapshot_id FROM phenotype_snapshot")] == [old]


@pytest.mark.parametrize("algorithm", ["identity", "king-robust"])
def test_genetic_parity_cli_and_cohort_oracle(bundle, algorithm):
    db, sidecar, manifest = bundle
    snapshot = import_snapshot(*bundle)["snapshot_id"]
    with closing(read_only(db)) as conn:
        pair = asdict(match_pair(conn, "S0", "S1", algorithm=algorithm))
        results = match_against_all(conn, "S0", top_n=2, algorithm=algorithm)
    report = match_report(db, sidecar, snapshot, "S0", "S1", algorithm=algorithm)
    assert report["genetic"] == pair
    report = match_report(db, sidecar, snapshot, "S0", algorithm=algorithm, top_n=2)
    expected = asdict(results) if algorithm == "king-robust" else [asdict(r) for r in results]
    assert report["genetic"] == expected
    predicate = {"ontology": "HPO", "term_id": "HP:0000001", "presence": "Absent"}
    report = match_report(
        db, sidecar, snapshot, "S0", algorithm=algorithm, top_n=1, predicate=predicate
    )
    assert report["cohort"]["sample_ids"] == ["S4"]
    assert report["cohort"]["query_outside_cohort"]
    with closing(read_only(db)) as conn:
        oracle = asdict(match_pair(conn, "S0", "S4", algorithm=algorithm))
    values = report["genetic"]["matches"] if algorithm == "king-robust" else report["genetic"]
    assert values == [oracle]
    result = CliRunner().invoke(
        cli,
        [
            "match-sample",
            "S0",
            "--db",
            str(db),
            "--algorithm",
            algorithm,
            "--phenotype-db",
            str(sidecar),
            "--phenotype-snapshot",
            snapshot,
            "--phenotype-term",
            "HPO",
            "HP:0000001",
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["cohort"]["count"] == 3


def test_no_sidecar_access_by_default(bundle, monkeypatch):
    db, _, _ = bundle
    monkeypatch.setattr(
        Phenotypes, "__init__", lambda *a, **k: pytest.fail("Unexpected clinical I/O")
    )
    for command in [("match-samples", "S0", "S1"), ("match-sample", "S0")]:
        result = CliRunner().invoke(cli, [*command, "--db", str(db), "--json"])
        assert result.exit_code == 0, result.output


def test_binding_detects_changed_genotypes_and_all_pairs_stream(bundle, tmp_path):
    db, sidecar, _ = bundle
    snapshot = import_snapshot(*bundle)["snapshot_id"]
    out = tmp_path / "export"
    report = export_pairs(db, sidecar, snapshot, out)
    assert report["pairs"] == 15
    assert len((out / "pairs.jsonl").read_text().splitlines()) == 15
    assert len((out / "participants.jsonl").read_text().splitlines()) == 4
    assert len((out / "observations.jsonl").read_text().splitlines()) == 6
    for name, sha in report["files"].items():
        assert file_sha256(out / name) == sha
    with sqlite3.connect(db) as conn:
        conn.execute('UPDATE genotype_call SET dosage=2 WHERE sample_id="S0"')
    with pytest.raises(ValueError, match="artifact"):
        Phenotypes(sidecar, snapshot, db)


def test_batched_reads_and_readonly(bundle):
    db, sidecar, _ = bundle
    snapshot = import_snapshot(*bundle)["snapshot_id"]
    with closing(Phenotypes(sidecar, snapshot, db)) as reader:
        trace = []
        reader.conn.set_trace_callback(trace.append)
        reader.enrich(["S0", "S1", "S2", "S3"])
        assert len([s for s in trace if s.startswith("SELECT")]) == 3
        with pytest.raises(sqlite3.OperationalError):
            reader.conn.execute("DELETE FROM phenotype_assertion")
        reader.conn.setlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER, 3)
        assert len(reader.enrich(["S0", "S1", "S2", "S3"])["observations"]) == 4


def test_fixture_checksums():
    for name, sha in json.loads((FIXTURE / "checksums.json").read_text()).items():
        assert file_sha256(FIXTURE / name) == sha
    for p in FIXTURE.iterdir():
        if p.is_file():
            assert "gs://" not in p.read_text() and "X-Goog-Signature" not in p.read_text()


def test_import_interruption_input_mutation_and_duplicate_counts(bundle, monkeypatch):
    import vrs_matcher.phenotypes as module

    db, sidecar, manifest = bundle
    old = import_snapshot(*bundle)["snapshot_id"]
    append(
        manifest.parent / "phenotype.tsv", "P0\tHP:0000001\tPresent\tHPO\tInvented test assertion\n"
    )
    refresh(manifest)
    real_rows = module.rows

    def interrupted(path, *args):
        for number, row in real_rows(path, *args):
            yield number, row
            if Path(path).name == "phenotype.tsv":
                raise KeyboardInterrupt

    monkeypatch.setattr(module, "rows", interrupted)
    with pytest.raises(KeyboardInterrupt):
        import_snapshot(*bundle)
    with closing(read_only(sidecar)) as conn:
        assert conn.execute("SELECT COUNT(*) FROM phenotype_snapshot").fetchone()[0] == 1
    monkeypatch.setattr(module, "rows", real_rows)
    assert import_snapshot(*bundle)["repeated_assertions"] == 1
    append(manifest.parent / "phenotype.tsv", "P0\tHP:0000002\tPresent\tHPO\tx\n")
    refresh(manifest)

    def mutate(path, *args):
        yield from real_rows(path, *args)
        if Path(path).name == "mapping.tsv":
            append(manifest.parent / "phenotype.tsv", "\n")

    monkeypatch.setattr(module, "rows", mutate)
    with pytest.raises(ValueError, match="changed"):
        import_snapshot(*bundle)
    with closing(Phenotypes(sidecar, old, db)) as reader:
        assert reader.cohort("HPO", "HP:0000002")["count"] == 0


def test_empty_cohort_candidate_restrictions_and_cli_errors(bundle):
    db, sidecar, _ = bundle
    snapshot = import_snapshot(*bundle)["snapshot_id"]
    for algorithm in ("identity", "king-robust"):
        report = match_report(
            db,
            sidecar,
            snapshot,
            "S0",
            algorithm=algorithm,
            predicate={"ontology": "HPO", "term_id": "HP:9999999"},
        )
        assert report["genetic"] == (
            [] if algorithm == "identity" else {"matches": [], "unscorable": []}
        )
        report = match_report(
            db, sidecar, snapshot, "S0", "S1", algorithm=algorithm, candidate_vrs_ids=frozenset()
        )
        with closing(read_only(db)) as conn:
            expected = asdict(
                match_pair(conn, "S0", "S1", algorithm=algorithm, candidate_vrs_ids=frozenset())
            )
        assert report["genetic"] == expected
    with pytest.raises(ValueError, match="built-in"):
        match_report(
            db,
            sidecar,
            snapshot,
            "S0",
            plugin_file="not-used.py",
            predicate={"ontology": "HPO", "term_id": "HP:0000001"},
        )
    result = CliRunner().invoke(
        cli, ["match-sample", "S0", "--db", str(db), "--phenotype-snapshot", snapshot]
    )
    assert result.exit_code != 0 and "both" in result.output


def test_acquisition_atomicity_and_replay(tmp_path, monkeypatch):
    from vrs_matcher.gregor_inputs import stage_inputs

    source = tmp_path / "download.tsv"
    source.write_text("synthetic\n")
    manifest = tmp_path / "acquisition.json"
    data = {
        "release": "synthetic",
        "workspace": "invented",
        "consent_group": "synthetic",
        "model_version": "1.12",
        "files": [{"name": "phenotype.tsv", "path": source.name, "sha256": file_sha256(source)}],
    }
    manifest.write_text(json.dumps(data))
    out = tmp_path / "ready"
    assert stage_inputs(manifest, out) == stage_inputs(manifest, out)
    source.write_text("truncated")
    with pytest.raises(ValueError, match="Checksum"):
        stage_inputs(manifest, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()
    source.write_text("synthetic\n")

    def fail_copy(*a):
        raise OSError("interrupted")

    monkeypatch.setattr(shutil, "copyfile", fail_copy)
    with pytest.raises(OSError):
        stage_inputs(manifest, tmp_path / "interrupted")
    assert not (tmp_path / "interrupted").exists()


def test_acquisition_requires_vcf_index(tmp_path):
    from vrs_matcher.gregor_inputs import stage_inputs

    source = FIXTURE / "calls.vcf"
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "release": "synthetic",
                "workspace": "test",
                "consent_group": "synthetic",
                "model_version": "1.11",
                "files": [
                    {
                        "name": "calls.vcf",
                        "path": str(source.resolve()),
                        "sha256": file_sha256(source),
                        "kind": "vcf",
                        "index": "calls.vcf.tbi",
                    }
                ],
            }
        )
    )
    with pytest.raises(ValueError, match="TBI/CSI"):
        stage_inputs(manifest, tmp_path / "output")
    assert not (tmp_path / "output").exists()
