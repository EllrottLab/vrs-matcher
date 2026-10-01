"""Generate wholly invented GREGoR inputs; never downloads or reads clinical data."""

import argparse
import csv
import json
from pathlib import Path

from vrs_matcher.genotypes import file_sha256


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def generate(output, participants=6, assertions=1, version="1.12"):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output / "participant.tsv",
        ["participant_id", "family_id", "affected_status"],
        (
            {"participant_id": f"P{i}", "family_id": f"F{i // 2}", "affected_status": "Unknown"}
            for i in range(participants)
        ),
    )
    fields = ["participant_id", "term_id", "presence", "ontology", "additional_details"]
    write_tsv(
        output / "phenotype.tsv",
        fields,
        (
            {
                "participant_id": f"P{i}",
                "term_id": f"HP:{j + 1:07}",
                "presence": ["Present", "Absent", "Unknown"][i % 3],
                "ontology": "HPO",
                "additional_details": "Invented test assertion",
            }
            for i in range(participants - 1)
            for j in range(assertions)
        ),
    )
    write_tsv(
        output / "mapping.tsv",
        ["sample_id", "participant_id"],
        (
            {"sample_id": f"S{i}", "participant_id": f"P{0 if i == 1 else i}"}
            for i in range(participants - 1)
        ),
    )
    source = {"source_id": "synthetic"}
    for table in ["participant", "phenotype", "mapping"]:
        path = output / f"{table}.tsv"
        source[table] = {"path": path.name, "sha256": file_sha256(path)}
    manifest = {
        "model_version": version,
        "genetic": {"artifact_id": "synthetic-king", "sha256": "SET_AFTER_INDEX_BUILD"},
        "sources": [source],
        "provenance": {
            "synthetic": True,
            "generator_version": 1,
            "participants": participants,
            "assertions_per_person": assertions,
        },
    }
    (output / "bundle.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Expand the existing hand-calculated KING fixture, preserving missing and reference calls.
    oracle = Path(__file__).resolve().parents[1] / "tests/data/king/oracle.vcf"
    with (output / "calls.vcf").open("w") as stream:
        for line in oracle.read_text().splitlines():
            if line.startswith("##"):
                stream.write(line + "\n")
            elif line.startswith("#CHROM"):
                stream.write(
                    "\t".join(line.split("\t")[:9] + [f"S{i}" for i in range(participants)]) + "\n"
                )
            else:
                fields = line.split("\t")
                stream.write(
                    "\t".join(fields[:9] + [fields[9 + i % 2] for i in range(participants)]) + "\n"
                )
    (output / "panel.tsv").write_bytes(oracle.with_name("panel.tsv").read_bytes())
    hashes = {
        p.name: file_sha256(p)
        for p in sorted(output.iterdir())
        if p.is_file()
        and p.name
        in [
            "calls.vcf",
            "panel.tsv",
            "participant.tsv",
            "phenotype.tsv",
            "mapping.tsv",
            "bundle.json",
        ]
    }
    (output / "checksums.json").write_text(json.dumps(hashes, indent=2) + "\n")
    return output / "bundle.json"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--participants", type=int, default=6)
    parser.add_argument("--assertions", type=int, default=1)
    parser.add_argument("--version", choices=["1.11", "1.12"], default="1.12")
    args = parser.parse_args()
    if args.participants < 2 or args.assertions < 1:
        parser.error("Require at least 2 participants and 1 assertion.")
    generate(args.output, args.participants, args.assertions, args.version)
