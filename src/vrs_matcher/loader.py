"""VCF ingestion utilities for building the sample-allele index.

This module parses VRS-annotated VCF records and transforms per-sample genotype
calls into normalized rows suitable for persistence in SQLite.
"""

from collections.abc import Iterator
from pathlib import Path

import cyvcf2

from .db import insert_alleles, open_db, register_samples
from .models import Zygosity

DEFAULT_GQ_THRESHOLD: float = 20
DEFAULT_DP_THRESHOLD: int = 0

# cyvcf2 uses INT_MIN for missing integer FORMAT fields.
_MISSING_INT = -2147483648


def _zygosity(gt_alleles: tuple[int | None, ...]) -> Zygosity:
    """Infer normalized zygosity from GT allele indexes.

    Args:
        gt_alleles: Allele index tuple from a sample genotype where ``0`` is
            reference, positive values map to ALT alleles, and ``None`` or
            negative values represent missing calls.

    Returns:
        The inferred :class:`vrs_matcher.models.Zygosity` value.
    """

    if all(a is None or a < 0 for a in gt_alleles):
        return Zygosity.NO_CALL
    valid = [a for a in gt_alleles if a is not None and a >= 0]
    non_ref = [a for a in valid if a != 0]
    if not non_ref:
        return Zygosity.REF
    has_missing = any(a is None or a < 0 for a in gt_alleles)
    if not has_missing and len(set(non_ref)) == 1 and len(non_ref) == len(valid):
        return Zygosity.HOM_ALT
    return Zygosity.HET


def genotype_to_vrs_ids(vrs_ids: list[str], gt_alleles: tuple[int | None, ...]) -> list[str]:
    """Map genotype allele indexes to VRS allele identifiers.

    Handles multi-allelic sites, phased or unphased genotypes, and missing
    calls. Allele index ``0`` (REF) and missing values (``None`` or negative)
    are ignored.

    Args:
        vrs_ids: Ordered list of VRS IDs corresponding to ALT alleles in the
            record.
        gt_alleles: Allele index tuple for one sample.

    Returns:
        List of carried VRS IDs, one entry per carried ALT allele.
    """
    carried = []
    for allele_index in gt_alleles:
        if allele_index is None or allele_index < 0 or allele_index == 0:
            continue
        idx = allele_index - 1
        if idx < len(vrs_ids):
            carried.append(vrs_ids[idx])
    return carried


def _gt_string(gt_alleles: tuple[int | None, ...], phased: bool) -> str:
    """Format GT allele indexes as a VCF-style genotype string.

    Args:
        gt_alleles: Allele index tuple for one sample.
        phased: Whether alleles are phased.

    Returns:
        Formatted genotype string such as ``"0/1"``, ``"1|1"``, or ``"./1"``.
    """

    sep = "|" if phased else "/"
    parts = [str(a) if (a is not None and a >= 0) else "." for a in gt_alleles]
    return sep.join(parts)


def _format_scalar(arr, sample_idx: int) -> float | int | None:
    """Extract one scalar FORMAT value for a sample.

    Args:
        arr: FORMAT matrix returned by ``cyvcf2`` for a field (for example,
            GQ or DP), or ``None`` when the field is absent.
        sample_idx: Zero-based sample index.

    Returns:
        Scalar value for the sample, or ``None`` when missing/NaN.
    """

    if arr is None:
        return None
    val = arr[sample_idx][0]
    if isinstance(val, float) and val != val:  # NaN
        return None
    if val == _MISSING_INT:
        return None
    return val


def _iter_rows(
    vcf_path: str | Path,
    *,
    source_dataset: str | None = None,
    gq_threshold: float = DEFAULT_GQ_THRESHOLD,
    dp_threshold: int = DEFAULT_DP_THRESHOLD,
    include_no_call: bool = False,
    candidate_vrs_ids: set[str] | None = None,
) -> Iterator[tuple]:
    """Generate normalized sample-allele rows from a VRS-annotated VCF.

    Args:
        vcf_path: Path to VCF/BCF input annotated with ``VRS_Allele_IDs``.
        source_dataset: Optional dataset label to store in emitted rows.
        gq_threshold: Minimum GQ to include when GQ is present.
        dp_threshold: Minimum DP to include when DP is present.
        include_no_call: Whether NO_CALL genotypes should be considered.
        candidate_vrs_ids: Optional allowlist of VRS IDs.

    Yields:
        Tuples in ``sample_allele`` insert order:
        ``(sample_id, vrs_id, gt, zygosity, chrom, pos, gq, dp, source_dataset)``.

    Raises:
        OSError: If the VCF input cannot be opened.
        Exception: Propagates parsing errors from ``cyvcf2``.
    """

    vcf = cyvcf2.VCF(str(vcf_path))
    sample_names: list[str] = vcf.samples

    try:
        for record in vcf:
            # Only include PASS (or unfiltered) records.
            f = record.FILTER
            if f and f != "PASS":
                continue

            vrs_ids_raw = record.INFO.get("VRS_Allele_IDs")
            if not vrs_ids_raw:
                continue
            vrs_ids = list(vrs_ids_raw)

            chrom = record.CHROM
            pos = record.POS
            gq_arr = record.format("GQ")
            dp_arr = record.format("DP")

            for sample_idx, sample_id in enumerate(sample_names):
                raw_gt = record.genotypes[sample_idx]  # [a1, a2, ..., phased_bool]
                phased = bool(raw_gt[-1])
                alleles: tuple[int | None, ...] = tuple(None if a < 0 else a for a in raw_gt[:-1])

                zyg = _zygosity(alleles)
                if zyg == Zygosity.NO_CALL and not include_no_call:
                    continue
                if zyg == Zygosity.REF:
                    continue

                sample_gq = _format_scalar(gq_arr, sample_idx)
                sample_dp = _format_scalar(dp_arr, sample_idx)

                if sample_gq is not None and sample_gq < gq_threshold:
                    continue
                if sample_dp is not None and sample_dp < dp_threshold:
                    continue

                gt_str = _gt_string(alleles, phased)
                # The sample-allele index stores carried alleles only. A
                # NO_CALL genotype does not identify any carried ALT allele, so
                # do not emit allele rows for it; otherwise downstream matching
                # logic can incorrectly treat the sample as carrying every ALT
                # VRS ID at this locus.
                if zyg == Zygosity.NO_CALL:
                    continue

                carried = genotype_to_vrs_ids(vrs_ids, alleles)

                for vrs_id in dict.fromkeys(carried):  # one row per (sample_id, vrs_id)
                    if candidate_vrs_ids is not None and vrs_id not in candidate_vrs_ids:
                        continue
                    yield (
                        sample_id,
                        vrs_id,
                        gt_str,
                        zyg.value,
                        chrom,
                        pos,
                        float(sample_gq) if sample_gq is not None else None,
                        int(sample_dp) if sample_dp is not None else None,
                        source_dataset,
                    )
    finally:
        vcf.close()


def load_samples(
    vcf_path: str | Path,
    db_path: str | Path,
    *,
    source_dataset: str | None = None,
    gq_threshold: float = DEFAULT_GQ_THRESHOLD,
    dp_threshold: int = DEFAULT_DP_THRESHOLD,
    include_no_call: bool = False,
    candidate_vrs_ids: set[str] | None = None,
) -> int:
    """Load a VRS-annotated VCF into the SQLite sample-allele index.

    Args:
        vcf_path: Path to VRS-annotated VCF input.
        db_path: Path to SQLite database file.
        source_dataset: Optional dataset label persisted with each row.
        gq_threshold: Minimum GQ threshold when GQ exists.
        dp_threshold: Minimum DP threshold when DP exists.
        include_no_call: Whether to include NO_CALL genotypes.
        candidate_vrs_ids: Optional allowlist of VRS IDs.

    Returns:
        Number of allele rows inserted.
    """

    _BATCH = 10_000
    conn = open_db(db_path)

    # Register every sample in the VCF header so that samples with zero
    # surviving alleles after filtering are still visible in the index.
    _vcf_header = cyvcf2.VCF(str(vcf_path))
    register_samples(conn, _vcf_header.samples)
    _vcf_header.close()

    total = 0
    batch: list[tuple] = []
    try:
        for row in _iter_rows(
            vcf_path,
            source_dataset=source_dataset,
            gq_threshold=gq_threshold,
            dp_threshold=dp_threshold,
            include_no_call=include_no_call,
            candidate_vrs_ids=candidate_vrs_ids,
        ):
            batch.append(row)
            if len(batch) >= _BATCH:
                insert_alleles(conn, batch)
                total += len(batch)
                batch = []
        if batch:
            insert_alleles(conn, batch)
            total += len(batch)
    finally:
        conn.close()
    return total
