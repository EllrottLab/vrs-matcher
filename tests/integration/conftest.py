"""Fixtures for real-data 1000 Genomes integration tests."""

import os
from pathlib import Path

import pytest

VCF_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
    "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
)

AFR_SAMPLES = [
    "NA18486",
    "NA18488",
    "NA18489",
    "NA18498",
    "NA18499",
    "NA18501",
    "NA18502",
    "NA18504",
    "NA18505",
    "NA18507",
]

EUR_SAMPLES = [
    "HG00096",
    "HG00097",
    "HG00099",
    "HG00100",
    "HG00101",
    "HG00102",
    "HG00103",
    "HG00105",
    "HG00106",
    "HG00107",
]

SAMPLE_TO_SUPER_POP = {**{sid: "AFR" for sid in AFR_SAMPLES}, **{sid: "EUR" for sid in EUR_SAMPLES}}


@pytest.fixture(scope="session")
def sample_to_super_population() -> dict[str, str]:
    """Return deterministic sample-to-super-population labels for assertions."""

    return SAMPLE_TO_SUPER_POP.copy()


@pytest.fixture(scope="session")
def annotated_vcf(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Create a local chr22 slice VCF with VRS allele IDs in INFO.

    Requires the environment variable ``GA4GH_VRS_DATAPROXY_URI`` to point at
    a reachable SeqRepo data source, for example::

        # local seqrepo install (recommended for CI)
        GA4GH_VRS_DATAPROXY_URI=seqrepo+file:///usr/local/share/seqrepo/2024-12-20

        # local seqrepo REST server
        GA4GH_VRS_DATAPROXY_URI=seqrepo+http://localhost:5000/seqrepo
    """

    pysam = pytest.importorskip("pysam", reason="integration tests need pysam for VCF slicing")
    pytest.importorskip(
        "ga4gh.vrs.extras.annotator.vcf",
        reason="integration tests require ga4gh.vrs[extras]",
    )

    dataproxy_uri = os.environ.get("GA4GH_VRS_DATAPROXY_URI")
    if not dataproxy_uri:
        pytest.skip(
            "Set GA4GH_VRS_DATAPROXY_URI to run integration tests, e.g. "
            "GA4GH_VRS_DATAPROXY_URI=seqrepo+file:///usr/local/share/seqrepo/2024-12-20"
        )

    from ga4gh.vrs.dataproxy import create_dataproxy
    from ga4gh.vrs.extras.annotator.vcf import VcfAnnotator

    try:
        proxy = create_dataproxy(dataproxy_uri)
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"Could not initialise SeqRepo data proxy ({dataproxy_uri!r}): {exc}")

    workdir = tmp_path_factory.mktemp("integration_1kg")
    subset_vcf = workdir / "chr22_subset.vcf.gz"
    annotated = workdir / "chr22_subset.annotated.vcf.gz"

    selected_samples = AFR_SAMPLES + EUR_SAMPLES

    try:
        with pysam.VariantFile(VCF_URL) as src:
            missing = [sid for sid in selected_samples if sid not in src.header.samples]
            if missing:
                pytest.skip(f"Selected 1KGP samples missing from upstream panel: {missing}")

            src.subset_samples(selected_samples)
            with pysam.VariantFile(str(subset_vcf), "wz", header=src.header) as out:
                for record in src.fetch("chr22", 16_000_000, 17_000_000):
                    out.write(record)
    except OSError as exc:
        pytest.skip(f"Unable to fetch 1KGP remote VCF slice: {exc}")

    annotator = VcfAnnotator(data_proxy=proxy)
    annotator.annotate(subset_vcf, annotated)

    return annotated

