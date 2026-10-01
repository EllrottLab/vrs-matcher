"""Deterministic offline phenotype phase benchmark; no biological performance claims."""

import argparse
import json
import platform
import resource
import statistics
import time
from contextlib import closing
from pathlib import Path

from generate_gregor_fixture import generate

from vrs_matcher.genotypes import file_sha256, load_genotypes
from vrs_matcher.matcher import match_against_all
from vrs_matcher.phenotypes import (
    Phenotypes,
    canonical,
    import_snapshot,
    read_only,
    snapshot_diff,
)


def run(output, participants, assertions, repeats):
    output.mkdir(parents=True, exist_ok=True)
    inputs = output / "inputs"
    manifest = generate(inputs, participants, assertions)
    db = output / "genetic.db"
    started = time.perf_counter()
    load_genotypes(inputs / "calls.vcf", db, inputs / "panel.tsv")
    index_seconds = time.perf_counter() - started
    data = json.loads(manifest.read_text())
    data["genetic"]["sha256"] = file_sha256(db)
    manifest.write_text(json.dumps(data))
    timings = []
    for repetition in range(repeats + 1):
        times = {}
        started = time.perf_counter()
        sidecar = output / f"clinical-{repetition}.db"
        snapshot = import_snapshot(db, sidecar, manifest)["snapshot_id"]
        times["import"] = time.perf_counter() - started
        with closing(Phenotypes(sidecar, snapshot, db)) as clinical, closing(read_only(db)) as conn:
            for algorithm in ["identity", "king-robust"]:
                started = time.perf_counter()
                match_against_all(conn, "S0", top_n=20, algorithm=algorithm)
                times[algorithm + "_genetic_query"] = time.perf_counter() - started
            started = time.perf_counter()
            clinical.cohort("HPO", "HP:0000001")
            times["cohort"] = time.perf_counter() - started
            trace = []
            clinical.conn.set_trace_callback(trace.append)
            started = time.perf_counter()
            report = clinical.enrich([f"S{i}" for i in range(min(20, participants))])
            times["enrich_20"] = time.perf_counter() - started
            clinical.conn.set_trace_callback(None)
            query_count = len(trace)
            started = time.perf_counter()
            payload = canonical(report)
            times["serialize_20"] = time.perf_counter() - started
            plan = [
                list(r)
                for r in clinical.conn.execute(
                    "EXPLAIN QUERY PLAN "
                    "SELECT participant_id FROM phenotype_assertion WHERE snapshot_id=? "
                    "AND ontology=? AND term_id=? AND presence=?",
                    (snapshot, "HPO", "HP:0000001", "Present"),
                )
            ]
        started = time.perf_counter()
        snapshot_diff(sidecar, snapshot, snapshot)
        times["unchanged_diff"] = time.perf_counter() - started
        timings.append(times)
    result = {
        "participants": participants,
        "assertions_per_person": assertions,
        "model_version": "1.12",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "first_run": timings[0],
        "warm": {
            k: {
                "median": statistics.median(t[k] for t in timings[1:]),
                "min": min(t[k] for t in timings[1:]),
                "max": max(t[k] for t in timings[1:]),
            }
            for k in timings[0]
        },
        "repeats": repeats,
        "index_seconds": index_seconds,
        "enrichment_sql_count": query_count,
        "cohort_query_plan": plan,
        "output_bytes": len(payload.encode()),
        "sidecar_bytes": sidecar.stat().st_size,
        "genetic_bytes": db.stat().st_size,
        "peak_rss_native_units": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "process_cpu_seconds": time.process_time(),
        "manifest_sha256": file_sha256(manifest),
        "cache": "OS cache uncontrolled; first run separate, subsequent runs warm",
        "limitation": "Eight-marker synthetic fixture; process RSS includes all phases.",
    }
    (output / "report.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("--participants", type=int, default=1000)
    p.add_argument("--assertions", type=int, default=10)
    p.add_argument("--repeats", type=int, default=5)
    args = p.parse_args()
    if args.output.exists() or args.participants < 2 or args.assertions < 1 or args.repeats < 1:
        p.error("Choose a new output directory, >=2 participants, >=1 assertion and >=1 repeat.")
    run(args.output, args.participants, args.assertions, args.repeats)
