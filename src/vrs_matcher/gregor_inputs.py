"""Publish verified, locally downloaded GREGoR inputs; no cloud credentials involved."""

import json
import shutil
import tempfile
from pathlib import Path

import cyvcf2

from .genotypes import file_sha256
from .phenotypes import canonical, contract


def stage_inputs(manifest_path, output):
    """Copy a manifest's files into a new atomic bundle, checking bytes and VCF indexes."""
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    contract(manifest["model_version"])
    if not all(manifest.get(k) for k in ("release", "workspace", "consent_group")):
        raise ValueError("Acquisition requires release, workspace and consent_group.")
    files = manifest["files"]
    names = [f["name"] for f in files]
    if (
        not names
        or len(names) != len(set(names))
        or any(Path(n).name != n or n in (".", "..", "manifest.json") for n in names)
    ):
        raise ValueError("Unique plain artifact filenames are required.")
    destination = Path(output)
    if destination.exists():
        existing = json.loads((destination / "manifest.json").read_text())
        expected = {k: v for k, v in manifest.items() if k != "files"}
        expected["files"] = [{k: v for k, v in f.items() if k != "path"} for f in files]
        if existing != expected or any(
            file_sha256(destination / f["name"]) != f["sha256"] for f in files
        ):
            raise ValueError("Existing acquisition differs; choose a new destination.")
        return existing
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        stage = Path(tmp) / "bundle"
        stage.mkdir()
        for spec in files:
            source = manifest_path.parent / spec["path"]
            target = stage / spec["name"]
            shutil.copyfile(source, target)
            if file_sha256(target) != spec["sha256"]:
                raise ValueError(f"Checksum mismatch: {spec['name']}")
        for spec in files:
            if spec.get("kind") == "vcf":
                index = spec.get("index")
                if index not in names or index not in (
                    spec["name"] + ".tbi",
                    spec["name"] + ".csi",
                ):
                    raise ValueError("VCF requires its declared matching TBI/CSI artifact.")
                vcf = cyvcf2.VCF(str(stage / spec["name"]))
                try:
                    if not vcf.samples or len(vcf.samples) != len(set(vcf.samples)):
                        raise ValueError("VCF requires unique sample genotypes, not sites only.")
                    if not vcf.seqnames:
                        raise ValueError("VCF must declare contigs.")
                    next(iter(vcf(f"{vcf.seqnames[0]}:1-1")), None)
                finally:
                    vcf.close()
        published = {k: v for k, v in manifest.items() if k != "files"}
        published["files"] = [{k: v for k, v in f.items() if k != "path"} for f in files]
        (stage / "manifest.json").write_text(canonical(published) + "\n")
        stage.rename(destination)
    return published
