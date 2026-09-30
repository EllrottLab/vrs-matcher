"""Opt-in, atomic indexing of called genotypes on a declared VRS SNP panel."""

import csv
import hashlib
import json
import math
import re
from collections import Counter
from itertools import batched
from pathlib import Path

import cyvcf2

from .db import insert_alleles, open_db, register_samples
from .loader import _format_scalar, _iter_rows

_VRS = re.compile(r"ga4gh:VA\.[A-Za-z0-9_-]{32}\Z")
_SEQUENCE = re.compile(r"ga4gh:SQ\.[A-Za-z0-9_-]{32}\Z")
_COLUMNS = {
    "reference",
    "annotation_version",
    "chrom",
    "sequence_id",
    "pos",
    "ref",
    "alt",
    "vrs_id",
}


def file_sha256(path: str | Path) -> str:
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_panel(path: str | Path) -> tuple[dict, dict, dict]:
    """Read a TSV with one explicit reference sequence identity per autosome."""
    markers, aliases = {}, {}
    loci = set()
    metadata = None
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if set(reader.fieldnames or []) != _COLUMNS or len(reader.fieldnames) != len(_COLUMNS):
            raise ValueError(f"Panel columns must be: {', '.join(sorted(_COLUMNS))}")
        for row in reader:
            if None in row or any(not v or v != v.strip() for v in row.values()):
                raise ValueError("Panel cells must be nonempty, without surrounding whitespace.")
            current = {k: row[k] for k in ("reference", "annotation_version")}
            if metadata is not None and current != metadata:
                raise ValueError("Panel reference and annotation_version must be consistent.")
            metadata = current
            chrom = row["chrom"].removeprefix("chr")
            if chrom not in {str(i) for i in range(1, 23)}:
                raise ValueError("Panel must contain human autosomes 1–22 only.")
            seq = row["sequence_id"]
            if not _SEQUENCE.fullmatch(seq) or not _VRS.fullmatch(row["vrs_id"]):
                raise ValueError(
                    "Panel requires ga4gh:SQ and ga4gh:VA identifiers with 32-char digests."
                )
            if (
                row["ref"] not in "ACGT"
                or row["alt"] not in "ACGT"
                or (len(row["ref"]) != 1 or len(row["alt"]) != 1 or row["ref"] == row["alt"])
            ):
                raise ValueError(
                    "Panel markers must be biallelic SNVs with distinct A/C/G/T alleles."
                )
            pos = int(row["pos"])
            if pos < 1 or (seq, pos) in loci or row["vrs_id"] in markers:
                raise ValueError("Invalid position or duplicate panel marker/locus.")
            for alias in (chrom, f"chr{chrom}"):
                if alias in aliases and aliases[alias] != seq:
                    raise ValueError("Conflicting sequence identity for panel chromosome.")
                aliases[alias] = seq
            if seq in aliases.values() and any(
                a.removeprefix("chr") != chrom and s == seq for a, s in aliases.items()
            ):
                raise ValueError("A reference sequence cannot identify two different autosomes.")
            markers[row["vrs_id"]] = (seq, pos, row["ref"], row["alt"])
            loci.add((seq, pos))
    if not markers:
        raise ValueError("Panel must contain at least one marker.")
    metadata["panel_id"] = file_sha256(path)
    return metadata, markers, aliases


def load_genotypes(vcf_path, db_path, panel_path, *, source_dataset=None, gq=20, dp=0) -> dict:
    """Load both indexes in one transaction; absent call rows always mean unknown.

    VCF headers must declare ##reference, ##vrs_annotation, and a refget
    attribute on each used contig. Digests are checked against the panel;
    sequence/allele normalization itself remains an upstream responsibility.
    """
    if not math.isfinite(gq) or not math.isfinite(dp) or gq < 0 or dp < 0:
        raise ValueError("GQ and DP thresholds must be finite and nonnegative.")
    metadata, markers, aliases = read_panel(panel_path)
    metadata["qc"] = {"gq": gq, "dp": dp, "missing_qc": "allow", "filter": "PASS_or_unfiltered"}
    metadata_json = json.dumps(metadata, sort_keys=True)
    provenance = {"input_sha256": file_sha256(vcf_path), "source_dataset": source_dataset}
    vcf = cyvcf2.VCF(str(vcf_path))
    conn = open_db(db_path)
    try:
        samples = vcf.samples
        if not samples or len(set(samples)) != len(samples):
            raise ValueError("VCF must contain unique, nonempty observation IDs.")
        for key, value in (
            ("reference", metadata["reference"]),
            ("vrs_annotation", metadata["annotation_version"]),
        ):
            values = [
                line.split("=", 1)[1]
                for line in vcf.raw_header.splitlines()
                if line.startswith(f"##{key}=")
            ]
            if values != [value]:
                raise ValueError(f"VCF ##{key} must match the panel: {value}")
        contigs = {}
        for header in vcf.header_iter():
            info = header.info()
            if info.get("HeaderType") == "CONTIG":
                try:
                    contigs[info["ID"]] = header["refget"].strip('"')
                except KeyError:
                    contigs[info["ID"]] = ""
        loci = {m[:2]: v for v, m in markers.items()}
        summary = Counter(observations=len(samples), panel_markers=len(markers))
        per_sample = {s: Counter() for s in samples}
        seen = set()
        calls = []
        # Reserve the writer before checking observation IDs, preventing append races.
        conn.execute("BEGIN IMMEDIATE")
        with conn:
            existing = conn.execute("SELECT metadata FROM genotype_panel").fetchone()
            if existing and json.loads(existing[0]) != metadata:
                raise ValueError(
                    "Incompatible panel/reference/annotation or QC policy; rebuild index."
                )
            registered = {r[0] for r in conn.execute("SELECT sample_id FROM samples")}
            if registered & set(samples):
                raise ValueError(
                    "Observation ID already exists; rename observations or rebuild index."
                )
            if not existing:
                conn.execute(
                    "INSERT INTO genotype_panel VALUES (1, ?, ?)",
                    (metadata["panel_id"], metadata_json),
                )
                conn.executemany(
                    "INSERT INTO genotype_marker VALUES (?, ?, ?, ?, ?)",
                    ((v, *m) for v, m in markers.items()),
                )
            register_samples(conn, samples, commit=False)
            for record in vcf:
                summary["records"] += 1
                seq = aliases.get(record.CHROM)
                if seq is None:
                    summary["excluded_nonpanel_contig_records"] += 1
                    continue
                if contigs.get(record.CHROM) != seq:
                    raise ValueError(f"Missing/inconsistent contig refget identity: {record.CHROM}")
                if (seq, record.POS) not in loci:
                    summary["excluded_off_panel_records"] += 1
                    continue
                vrs_id = loci[(seq, record.POS)]
                if vrs_id in seen:
                    raise ValueError(f"Duplicate panel record: {record.CHROM}:{record.POS}")
                seen.add(vrs_id)
                if (
                    len(record.ALT) != 1
                    or record.REF not in {"A", "C", "G", "T"}
                    or record.ALT[0] not in {"A", "C", "G", "T"}
                ):
                    summary["excluded_unsupported_records"] += 1
                    continue
                if (seq, record.POS, record.REF, record.ALT[0]) != markers[vrs_id]:
                    raise ValueError(f"Panel marker REF/ALT mismatch: {vrs_id}")
                if record.INFO.get("VRS_Allele_IDs") != vrs_id:
                    raise ValueError(
                        f"Missing/invalid VRS annotation or panel marker mismatch: {vrs_id}"
                    )
                if record.FILTER and record.FILTER != "PASS":
                    summary["excluded_filter_records"] += 1
                    continue
                summary["retained_records"] += 1
                formats = {}
                for key in ("GQ", "DP"):
                    try:
                        formats[key] = record.format(key)
                    except KeyError:
                        formats[key] = None
                genotypes = record.genotypes
                for i, sample in enumerate(samples):
                    stats = per_sample[sample]
                    stats["eligible_calls"] += 1
                    alleles = genotypes[i][:-1]
                    gq_value = _format_scalar(formats["GQ"], i)
                    dp_value = _format_scalar(formats["DP"], i)
                    for key, value in (("gq", gq_value), ("dp", dp_value)):
                        if value is None:
                            stats[f"missing_{key}"] += 1
                    if len(alleles) != 2:
                        reason = "non_diploid"
                    elif any(a < 0 for a in alleles):
                        reason = "missing_gt" if all(a < 0 for a in alleles) else "partial_gt"
                    elif any(a not in (0, 1) for a in alleles):
                        reason = "invalid_allele"
                    elif gq_value is not None and gq_value < gq:
                        reason = "low_gq"
                    elif dp_value is not None and dp_value < dp:
                        reason = "low_dp"
                    else:
                        dosage = sum(alleles)
                        stats[f"dosage_{dosage}"] += 1
                        calls.append((sample, vrs_id, dosage))
                        if len(calls) >= 10_000:
                            conn.executemany("INSERT INTO genotype_call VALUES (?, ?, ?)", calls)
                            calls.clear()
                        continue
                    stats[f"excluded_{reason}"] += 1
            conn.executemany("INSERT INTO genotype_call VALUES (?, ?, ?)", calls)
            summary["absent_panel_records"] = len(markers) - len(seen)
            # Reuse the identity projection unchanged, inside the same transaction.
            for batch in batched(
                _iter_rows(
                    vcf_path, source_dataset=source_dataset, gq_threshold=gq, dp_threshold=dp
                ),
                10_000,
                strict=False,
            ):
                insert_alleles(conn, batch, commit=False)
                summary["allele_rows"] += len(batch)
            if file_sha256(vcf_path) != provenance["input_sha256"]:
                raise ValueError("VCF changed during indexing; no observations were published.")
            for sample, stats in per_sample.items():
                conn.execute(
                    "INSERT INTO genotype_observation VALUES (?, 1, ?, ?)",
                    (
                        sample,
                        json.dumps(provenance),
                        json.dumps({"records": dict(summary), "calls": dict(stats)}),
                    ),
                )
            for stats in per_sample.values():
                summary.update(stats)
        return {"panel_id": metadata["panel_id"], **summary}
    finally:
        vcf.close()
        conn.close()
