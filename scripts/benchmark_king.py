#!/usr/bin/env python3
"""Opt-in, offline KING benchmark; fresh processes, raw artifacts, no extra dependencies."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import platform
import random
import shutil
import sqlite3
import statistics
import subprocess
import sys
import time
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path


def read_tsv(path):
    with open(path) as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows):
    rows = iter(rows)
    first = next(rows, None)
    if first is None:
        path.write_text("")
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(first), delimiter="\t")
        writer.writeheader()
        writer.writerow(first)
        writer.writerows(rows)


def digest(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def validate_prepared_input(args):
    """Refuse unfair comparisons: baseline input must already follow the panel/QC contract."""
    import cyvcf2

    from vrs_matcher.genotypes import read_panel
    from vrs_matcher.loader import _format_scalar

    _, markers, aliases = read_panel(args.panel)
    vcf = cyvcf2.VCF(str(args.vcf))
    seen = set()
    sample_count = len(vcf.samples)
    try:
        if len(vcf.samples) < 2:
            raise ValueError("Benchmark requires at least two observations.")
        for record in vcf:
            marker = record.INFO.get("VRS_Allele_IDs")
            if marker not in markers or marker in seen or len(record.ALT) != 1:
                raise ValueError("Benchmark VCF must contain unique markers from the panel only.")
            if markers[marker] != (
                aliases.get(record.CHROM),
                record.POS,
                record.REF,
                record.ALT[0],
            ):
                raise ValueError("Benchmark marker does not match panel.")
            seen.add(marker)
            if record.FILTER and record.FILTER != "PASS":
                raise ValueError("Prefilter failed records before benchmarking.")
            for key, threshold in [("GQ", 20), ("DP", 0)]:
                try:
                    arr = record.format(key)
                except KeyError:
                    arr = None
                for i in range(len(vcf.samples)):
                    val = _format_scalar(arr, i)
                    if val is not None and val < threshold:
                        raise ValueError(
                            "Materialize QC masking and remove low GQ/DP before benchmarking."
                        )
            for gt in record.genotypes:
                alleles = gt[:-1]
                if len(alleles) != 2 or not (
                    all(a in (0, 1) for a in alleles) or all(a == -1 for a in alleles)
                ):
                    raise ValueError("Benchmark permits complete diploid calls or ./. only.")
        if not seen:
            raise ValueError("Benchmark VCF has no panel records.")
    finally:
        vcf.close()
    return {"samples": sample_count, "markers": len(seen)}


def generate_synthetic_inputs(vcf_path, panel_path, *, samples, markers, seed):
    """Write reproducible synthetic diploid calls and a matching panel."""
    if samples < 2 or markers < 1:
        raise ValueError("Synthetic benchmark needs at least two samples and one marker.")
    rng = random.Random(seed)
    sample_ids = [f"S{i:05d}" for i in range(samples)]
    with panel_path.open("w") as panel, vcf_path.open("w") as vcf:
        panel.write("reference\tannotation_version\tchrom\tsequence_id\tpos\tref\talt\tvrs_id\n")
        vcf.write(
            "##fileformat=VCFv4.2\n"
            "##reference=synthetic-benchmark-reference\n"
            "##vrs_annotation=synthetic-benchmark-1\n"
            f"##contig=<ID=chr1,length={markers + 1},"
            "refget=ga4gh:SQ.AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA>\n"
            "##INFO=<ID=VRS_Allele_IDs,Number=A,Type=String,"
            'Description="Synthetic benchmark marker">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_ids) + "\n"
        )
        for marker in range(1, markers + 1):
            vrs_id = f"ga4gh:VA.{marker:032d}"
            panel.write(
                f"synthetic-benchmark-reference\tsynthetic-benchmark-1\t1\t"
                f"ga4gh:SQ.{'A' * 32}\t{marker}\tA\tC\t{vrs_id}\n"
            )
            genotypes = (("0/0", "0/1", "1/1")[rng.randrange(3)] for _ in sample_ids)
            vcf.write(
                f"chr1\t{marker}\t.\tA\tC\t.\tPASS\tVRS_Allele_IDs={vrs_id}\tGT:GQ:DP\t"
                + "\t".join(f"{gt}:30:20" for gt in genotypes)
                + "\n"
            )


def profile_indexed_pairs(context):
    """Time tiled SQLite retrieval separately from scoring/iteration."""
    from vrs_matcher.king import KingRobustPlugin

    class RetrievalTimer:
        def __init__(self, wrapped):
            self.wrapped = wrapped
            self.retrieval_seconds = 0.0

        def get_called_genotypes_many(self, *args, **kwargs):
            start = time.perf_counter()
            try:
                return self.wrapped.get_called_genotypes_many(*args, **kwargs)
            finally:
                self.retrieval_seconds += time.perf_counter() - start

        def __getattr__(self, name):
            return getattr(self.wrapped, name)

    timed_context = RetrievalTimer(context)
    start = time.perf_counter()
    pairs = scorable = checksum = 0
    for result in KingRobustPlugin().iter_all_pairs(timed_context):
        pairs += 1
        scorable += result.kinship is not None
        checksum += result.n_common + result.het_a + result.het_b + result.opposite_hom
    total_seconds = time.perf_counter() - start
    return {
        "retrieval_seconds": timed_context.retrieval_seconds,
        "scoring_iteration_seconds": max(0.0, total_seconds - timed_context.retrieval_seconds),
        "profile_seconds": total_seconds,
        "pairs": pairs,
        "scorable_pairs": scorable,
        "checksum": checksum,
    }


def worker(args):
    from vrs_matcher.db import list_samples, open_db
    from vrs_matcher.king import KingRobustPlugin
    from vrs_matcher.loader import load_samples
    from vrs_matcher.matcher import match_against_all
    from vrs_matcher.plugins import PluginContext

    if args.stage in ("load", "identity-load"):
        start = time.perf_counter()
        result = load_samples(
            args.vcf,
            args.db,
            index_genotypes=args.stage == "load",
            panel=args.panel if args.stage == "load" else None,
        )
        print(json.dumps({"index_build_seconds": time.perf_counter() - start, "summary": result}))
        return
    conn = open_db(args.db)
    try:
        samples = list_samples(conn)
        if args.stage == "query":
            print(json.dumps(asdict(match_against_all(conn, samples[0], algorithm="king-robust"))))
        elif args.stage == "profile":
            print(json.dumps(profile_indexed_pairs(PluginContext(conn))))
        else:
            pair_results = KingRobustPlugin().iter_all_pairs(PluginContext(conn))
            write_tsv(args.output / "pairs.tsv", (asdict(row) for row in pair_results))
    finally:
        conn.close()


def measure(command, log):
    """wait4 returns per-process resource usage, avoiding cumulative RSS across runs."""
    start = time.perf_counter()
    with log.open("w") as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
    return {
        "wall_seconds": time.perf_counter() - start,
        "cpu_seconds": usage.ru_utime + usage.ru_stime,
        "peak_rss_bytes": usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024),
        "exit_code": process.returncode,
    }


def summarize_workloads(timings):
    """Summarize index build, indexed matching, and direct VCFtools workloads."""
    repeats = {}
    for row in timings:
        if row["repeat"] == 0:
            continue
        repeats.setdefault(row["repeat"], {})[row["stage"]] = row["wall_seconds"]
    indexed_runs = {repeat: stages["all"] for repeat, stages in repeats.items() if "all" in stages}
    pipeline_runs = {
        repeat: stages["load"] + stages["all"]
        for repeat, stages in repeats.items()
        if "load" in stages and "all" in stages
    }
    vcftools_runs = {
        repeat: stages["vcftools"] for repeat, stages in repeats.items() if "vcftools" in stages
    }
    first_runs = {}
    for row in timings:
        if row["repeat"] == 0:
            first_runs[row["stage"]] = row["wall_seconds"]
    matched = sorted(indexed_runs.keys() & pipeline_runs.keys() & vcftools_runs.keys())

    def comparison(stage_runs, first_time):
        repeats_to_use = matched or ([0] if first_time is not None else [])
        if not repeats_to_use:
            return None
        if matched:
            plugin_times = [stage_runs[repeat] for repeat in matched]
            baseline_times = [vcftools_runs[repeat] for repeat in matched]
            basis = "subsequent runs"
        else:
            plugin_times = [first_time]
            baseline_times = [first_runs["vcftools"]]
            basis = "first run only"
        plugin_median = statistics.median(plugin_times)
        baseline_median = statistics.median(baseline_times)
        return {
            "plugin_median_seconds": plugin_median,
            "vcftools_median_seconds": baseline_median,
            "relative_performance_percent": (
                baseline_median / plugin_median * 100 if plugin_median else math.inf
            ),
            "basis": basis,
        }

    indexed_compare = comparison(indexed_runs, first_runs.get("all"))
    pipeline_compare = comparison(
        pipeline_runs,
        first_runs.get("load", 0) + first_runs.get("all", 0)
        if {"load", "all"} <= first_runs.keys()
        else None,
    )
    return {
        "indexed_all_pairs": indexed_compare,
        "index_plus_all_pairs": pipeline_compare,
        "index_build_median_seconds": statistics.median(
            [r["load"] for r in repeats.values() if "load" in r]
        )
        if any("load" in r for r in repeats.values())
        else first_runs.get("load"),
        "query_median_seconds": statistics.median(
            [r["query"] for r in repeats.values() if "query" in r]
        )
        if any("query" in r for r in repeats.values())
        else first_runs.get("query"),
    }


def summarize_labels(args, scores):
    """Summarize externally supplied labels, without deriving truth from either tool."""
    from vrs_matcher.genotypes import read_panel

    panel_size = len(read_panel(args.panel)[1])
    metrics = ["kinship", "kinship_within_family"]
    if all("vcftools_phi" in row for row in scores):
        metrics.append("vcftools_phi")
    samples = read_tsv(args.samples)
    if any(
        not {"sample", "donor", "family", "ancestry"} <= row.keys()
        or any(not row[k] for k in ("sample", "donor", "family", "ancestry"))
        for row in samples
    ):
        raise ValueError(
            "Sample labels require sample, donor, family, ancestry columns and values."
        )
    sample_map = {r["sample"]: r for r in samples}
    if len(sample_map) != len(samples):
        raise ValueError("Duplicate sample labels.")
    pairs = {tuple(sorted((r["sample_a"], r["sample_b"]))): r for r in scores}
    if set(sample_map) != {s for p in pairs for s in p}:
        raise ValueError("Sample labels must exactly match the benchmark cohort.")
    labels = read_tsv(args.relationships)
    required = {"duplicate", "parent_offspring", "siblings", "unrelated"}
    if not required <= {r["relationship"] for r in labels}:
        raise ValueError(
            "Biological evaluation requires duplicate, parent_offspring, siblings, unrelated."
        )
    groups, annotated = {}, []
    seen = set()
    for label in labels:
        key = tuple(sorted((label["sample_a"], label["sample_b"])))
        if key not in pairs or key in seen:
            raise ValueError(f"Unknown or duplicated labeled pair: {key}")
        seen.add(key)
        row = dict(pairs[key])
        a, b = (sample_map[s] for s in key)
        relation = label["relationship"]
        if relation == "duplicate" and a["donor"] != b["donor"]:
            raise ValueError(f"Duplicate label must share a donor: {key}")
        ancestry = "/".join(sorted({a["ancestry"], b["ancestry"]}))
        fraction = int(row["n_common"]) / panel_size
        overlap = ">=99%" if fraction >= 0.99 else "90–99%" if fraction >= 0.9 else "<90%"
        row.update(relationship=relation, ancestry=ancestry, common_fraction=fraction)
        annotated.append(row)
        groups.setdefault((relation, ancestry, overlap), []).append(row)
    write_tsv(args.output / "labeled-pairs.tsv", annotated)
    summaries = []
    for (relation, ancestry, overlap), rows in groups.items():
        for metric in metrics:
            values = [float(r[metric]) for r in rows if r[metric] != ""]
            summaries.append(
                {
                    "relationship": relation,
                    "ancestry": ancestry,
                    "callable_overlap": overlap,
                    "metric": metric,
                    "pairs": len(rows),
                    "scorable": len(values),
                    "mean": statistics.mean(values) if values else "",
                    "min_common_loci": min(int(r["n_common"]) for r in rows),
                    "max_common_loci": max(int(r["n_common"]) for r in rows),
                }
            )
    write_tsv(args.output / "relationship-summary.tsv", summaries)
    # Query-level duplicate retrieval. Count tied best hits separately from unique recovery.
    retrieval = []
    for sample, meta in sample_map.items():
        truth = {s for s, m in sample_map.items() if m["donor"] == meta["donor"] and s != sample}
        if not truth:
            continue
        for metric in metrics:
            ranked = []
            for key, row in pairs.items():
                if sample in key and row[metric] != "":
                    ranked.append((float(row[metric]), key[1] if key[0] == sample else key[0]))
            ranked.sort(key=lambda x: (-x[0], x[1]))
            best = {peer for score, peer in ranked if score == ranked[0][0]} if ranked else set()
            retrieval.append(
                {
                    "sample": sample,
                    "family": meta["family"],
                    "metric": metric,
                    "unique_correct": int(len(best) == 1 and bool(best & truth)),
                    "tied_best": int(len(best) > 1),
                    "unscorable": int(not ranked),
                }
            )
    write_tsv(args.output / "duplicate-retrieval.tsv", retrieval)
    # Family-block bootstrap of query accuracy; keep correlated queries together.
    rng = random.Random(0)
    intervals = {}
    for metric in metrics:
        blocks = {}
        for row in retrieval:
            if row["metric"] == metric:
                blocks.setdefault(row["family"], []).append(row["unique_correct"])
        if len(blocks) < 2:
            intervals[metric] = {"reason": "Fewer than two independent family blocks"}
            continue
        values = list(blocks.values())
        boot = sorted(
            statistics.mean(list(itertools.chain.from_iterable(rng.choices(values, k=len(values)))))
            for _ in range(1000)
        )
        intervals[metric] = {
            "family_blocks": len(blocks),
            "bootstrap_seed": 0,
            "accuracy_95_percent_interval": [boot[24], boot[974]],
        }
    (args.output / "uncertainty.json").write_text(json.dumps(intervals, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, help="Prefiltered, VRS-annotated VCF(.gz).")
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--output", type=Path, required=True, help="New output directory.")
    parser.add_argument("--vcftools", default="vcftools")
    parser.add_argument("--samples", type=Path, help="sample/donor/family/ancestry TSV")
    parser.add_argument("--relationships", type=Path, help="sample_a/sample_b/relationship TSV")
    parser.add_argument(
        "--smoke", action="store_true", help="Skip biological labels; no accuracy claim."
    )
    parser.add_argument(
        "--repeats", type=int, default=6, help="First run plus five repeats by default."
    )
    parser.add_argument(
        "--synthetic-samples",
        type=int,
        help="Generate deterministic synthetic input with this many observations.",
    )
    parser.add_argument(
        "--synthetic-markers",
        type=int,
        help="Generate deterministic synthetic input with this many SNP markers.",
    )
    parser.add_argument("--seed", type=int, default=0, help="Seed for synthetic input generation.")
    parser.add_argument(
        "--stage",
        choices=["load", "identity-load", "query", "all", "profile"],
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--db", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.stage:
        worker(args)
        return
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if not args.smoke and (not args.samples or not args.relationships):
        parser.error("Supply biological labels or explicitly select --smoke")
    synthetic_requested = args.synthetic_samples is not None or args.synthetic_markers is not None
    if synthetic_requested:
        if args.synthetic_samples is None or args.synthetic_markers is None:
            parser.error("Specify both --synthetic-samples and --synthetic-markers.")
        if args.vcf is not None or args.panel is not None:
            parser.error("Do not combine synthetic dimensions with --vcf or --panel.")
        if args.synthetic_samples < 2 or args.synthetic_markers < 1:
            parser.error("Synthetic input needs at least two samples and one marker.")
    elif args.vcf is None or args.panel is None:
        parser.error("Supply both --vcf and --panel, or use both synthetic dimensions.")
    executable = shutil.which(args.vcftools)
    if not executable:
        parser.error("VCFtools executable not found")
    args.output.mkdir(parents=True, exist_ok=False)
    synthetic_config = None
    if synthetic_requested:
        args.vcf = args.output / "synthetic.vcf"
        args.panel = args.output / "synthetic-panel.tsv"
        try:
            generate_synthetic_inputs(
                args.vcf,
                args.panel,
                samples=args.synthetic_samples,
                markers=args.synthetic_markers,
                seed=args.seed,
            )
        except ValueError as exc:
            parser.error(str(exc))
        synthetic_config = {
            "samples": args.synthetic_samples,
            "markers": args.synthetic_markers,
            "seed": args.seed,
        }
    dimensions = validate_prepared_input(args)
    root = Path(__file__).resolve().parents[1]
    manifest = {
        "command": sys.argv,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "processor": platform.processor(),
        "sqlite": sqlite3.sqlite_version,
        "python": sys.version,
        "mode": "smoke" if args.smoke else "labeled",
        "cohort": dimensions,
        "synthetic_input": synthetic_config,
        "cache": "uncontrolled OS cache; fresh processes; first run reported separately",
        "sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "dirty": subprocess.check_output(["git", "status", "--short"], cwd=root, text=True),
        "inputs": {str(p): digest(p) for p in [args.vcf, args.panel, root / "uv.lock"]},
        "packages": {p: version(p) for p in ["cyvcf2", "click"]},
        "vcftools": subprocess.check_output([executable, "--version"], text=True).strip(),
        "vcftools_sha256": digest(executable),
        "source_hashes": {
            str(p.relative_to(root)): digest(p)
            for p in [*sorted((root / "src").rglob("*.py")), Path(__file__).resolve()]
        },
        "preparation": "Inputs are prefiltered and annotated; upstream time is not measured.",
    }
    for path in (args.samples, args.relationships):
        if path:
            manifest["inputs"][str(path)] = digest(path)
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    timings = []
    for repeat in range(args.repeats):
        db = args.output / f"run-{repeat}.db"

        def command(stage, db=db):
            return [
                sys.executable,
                str(Path(__file__).resolve()),
                "--stage",
                stage,
                "--vcf",
                str(args.vcf),
                "--panel",
                str(args.panel),
                "--db",
                str(db),
                "--output",
                str(args.output),
            ]

        baseline = [
            executable,
            "--gzvcf" if args.vcf.suffix == ".gz" else "--vcf",
            str(args.vcf),
            "--relatedness2",
            "--out",
            str(args.output / f"vcftools-{repeat}"),
        ]
        stages = [
            ("load", command("load")),
            ("profile", command("profile")),
            ("query", command("query")),
            ("all", command("all")),
        ]
        stages.insert(len(stages) if repeat % 2 else 0, ("vcftools", baseline))
        for stage, cmd in stages:
            metrics = measure(cmd, args.output / f"{repeat}-{stage}.log")
            row = {
                "repeat": repeat,
                "stage": stage,
                **metrics,
                "index_build_seconds": "",
                "retrieval_seconds": "",
                "scoring_iteration_seconds": "",
                "profile_seconds": "",
                "samples": dimensions["samples"],
                "markers": dimensions["markers"],
                "pairs": "",
            }
            if stage == "load":
                load_result = json.loads((args.output / f"{repeat}-{stage}.log").read_text())
                row["index_build_seconds"] = load_result["index_build_seconds"]
            elif stage == "profile":
                phases = json.loads((args.output / f"{repeat}-{stage}.log").read_text())
                if phases["pairs"] != dimensions["samples"] * (dimensions["samples"] - 1) // 2:
                    raise RuntimeError("Profiled pair count does not match cohort dimensions.")
                row.update(
                    retrieval_seconds=phases["retrieval_seconds"],
                    scoring_iteration_seconds=phases["scoring_iteration_seconds"],
                    profile_seconds=phases["profile_seconds"],
                    pairs=phases["pairs"],
                )
            timings.append(row)
            write_tsv(args.output / "timings.tsv", timings)
            if metrics["exit_code"]:
                raise RuntimeError(f"{stage} failed; see {repeat}-{stage}.log and timings.tsv")
    # One separate identity-only build measures storage overhead, not speed superiority.
    args.db = args.output / "identity.db"
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--stage",
        "identity-load",
        "--vcf",
        str(args.vcf),
        "--panel",
        str(args.panel),
        "--db",
        str(args.db),
        "--output",
        str(args.output),
    ]
    metric = measure(cmd, args.output / "identity-load.log")
    if metric["exit_code"]:
        raise RuntimeError("Identity-only load failed; see identity-load.log")
    scores = read_tsv(args.output / "pairs.tsv")
    baseline_rows = {
        tuple(sorted((r["INDV1"], r["INDV2"]))): r
        for r in read_tsv(args.output / "vcftools-0.relatedness2")
    }
    comparison = []
    for row in scores:
        baseline_row = baseline_rows[tuple(sorted((row["sample_a"], row["sample_b"])))]
        phi = float(baseline_row["RELATEDNESS_PHI"])
        within = float(row["kinship_within_family"]) if row["kinship_within_family"] else None
        comparison.append(
            {
                **row,
                "vcftools_phi": phi if math.isfinite(phi) else "",
                "within_family_delta": within - phi
                if within is not None and math.isfinite(phi)
                else "",
                "vcftools_het_both": baseline_row["N_AaAa"],
                "vcftools_opposite_hom": baseline_row["N_AAaa"],
            }
        )
    write_tsv(args.output / "comparison.tsv", comparison)
    if not args.smoke:
        summarize_labels(args, comparison)
    report = [
        "# KING benchmark report",
        "",
        f"Mode: {manifest['mode']}. Cohort: {dimensions['samples']} observations, "
        f"{dimensions['markers']} markers, {len(scores)} unique pairs.",
        "",
        "Input preparation and VRS annotation are excluded; this is not end-to-end timing.",
        "VCFtools computes both pair directions and diagonals; the plugin all-pairs stage",
        "computes each unordered distinct pair once. Query is a separate workload.",
        "",
        "Workloads are reported separately: `load` builds the genotype and allele",
        "SQLite indexes; `query` matches one sample against the cohort; `all` matches",
        "all unordered distinct pairs using an existing index; `vcftools` reads the",
        "prepared VCF and computes `--relatedness2`. First-run time is repeat 0.",
        "Subsequent median and range exclude repeat 0. Peak RSS is the maximum process",
        "peak across all repeats for that stage, in MiB (1024² bytes). All times",
        "include process startup.",
        "",
        "| Stage | First run (s) | Subsequent median (s) | "
        "Subsequent range (s) | Max peak RSS (MiB) |",
        "|---|---:|---:|---|---:|",
    ]
    for stage in ["load", "query", "all", "vcftools"]:
        rows = [r for r in timings if r["stage"] == stage]
        walls = [r["wall_seconds"] for r in rows[1:]]
        median = f"{statistics.median(walls):.6f}" if walls else "NA"
        spread = f"{min(walls):.6f}–{max(walls):.6f}" if walls else "NA"
        report.append(
            f"| {stage} | {rows[0]['wall_seconds']:.6f} | {median} | "
            f"{spread} | {max(r['peak_rss_bytes'] for r in rows) / (1024**2):.2f} |"
        )
    report += [
        "",
        "## Isolated index, retrieval, and scoring phases",
        "",
        "These phase durations are measured inside the worker, excluding subprocess "
        "startup. `index build` times input parsing and database creation. `call "
        "retrieval` times SQLite reads into tiled genotype maps; `scoring + iteration` "
        "covers pair enumeration, score calculation, and result construction but "
        "excludes those retrieval calls. Process peak RSS is reported for the worker "
        "that measured each phase.",
        "",
        "| Phase | First run (s) | Subsequent median (s) | Subsequent range (s) | "
        "Max worker peak RSS (MiB) |",
        "|---|---:|---:|---|---:|",
    ]
    for label, stage, key in [
        ("Index build", "load", "index_build_seconds"),
        ("Call retrieval", "profile", "retrieval_seconds"),
        ("Scoring + pair iteration", "profile", "scoring_iteration_seconds"),
    ]:
        rows = [r for r in timings if r["stage"] == stage]
        values = [float(r[key]) for r in rows[1:]]
        median = f"{statistics.median(values):.6f}" if values else "NA"
        spread = f"{min(values):.6f}–{max(values):.6f}" if values else "NA"
        report.append(
            f"| {label} | {float(rows[0][key]):.6f} | {median} | {spread} | "
            f"{max(r['peak_rss_bytes'] for r in rows) / (1024**2):.2f} |"
        )
    workloads = summarize_workloads(timings)
    report += ["", "## Workload comparisons", ""]
    report.append(
        "Relative performance is `100 × VCFtools median / VRS-Matcher median`; "
        "VCFtools is 100%. Each VRS-Matcher median is computed from paired repeats."
    )
    report += [
        "",
        "| VRS-Matcher workload | VRS-Matcher median (s) | VCFtools direct median (s) | "
        "Relative performance (% of VCFtools speed) |",
        "|---|---:|---:|---:|",
    ]
    for label, summary in [
        ("Indexed all-pairs matching only", workloads["indexed_all_pairs"]),
        ("Index build + all-pairs matching", workloads["index_plus_all_pairs"]),
    ]:
        if summary is None:
            report.append(f"| {label} | NA | NA | NA |")
        else:
            report.append(
                f"| {label} | {summary['plugin_median_seconds']:.6f} | "
                f"{summary['vcftools_median_seconds']:.6f} | "
                f"{summary['relative_performance_percent']:.1f} |"
            )
    report += [
        "",
        "The first comparison measures matching against an already-built VRS-Matcher "
        "index; the VCFtools figure still includes reading the VCF, so it is not an "
        "algorithm-only comparison. The second includes VRS-Matcher index construction "
        "and is the one-shot pipeline view. VCFtools reads the prepared VCF directly. "
        "The tools also emit different pair sets as noted above. The `query` row is a "
        "separate one-versus-all workload and has no equivalent single-query VCFtools "
        "operation in this benchmark.",
    ]
    report += [
        "",
        f"Identity-only DB bytes: {args.db.stat().st_size}.",
        f"Genotype + allele DB bytes: {db.stat().st_size}.",
        f"Annotated VCF bytes: {args.vcf.stat().st_size}.",
        "",
        "See manifest.json, raw logs, pairs.tsv, and timings.tsv for reproducibility.",
        "No biological accuracy claim is made by a smoke run.",
    ]
    (args.output / "report.md").write_text("\n".join(report) + "\n")
    print(args.output / "report.md")


if __name__ == "__main__":
    main()
