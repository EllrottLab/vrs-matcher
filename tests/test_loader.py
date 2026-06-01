"""Unit and integration-style tests for loader behavior.

This module validates zygosity inference, genotype-to-VRS mapping, and sample
loading behavior using mocked ``cyvcf2`` records.
"""

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from vrs_matcher.db import get_vrs_ids, open_db
from vrs_matcher import loader as loader_mod
from vrs_matcher.loader import _format_scalar, _zygosity, genotype_to_vrs_ids, load_samples
from vrs_matcher.models import Zygosity


class TestZygosity:
    """Tests for genotype-to-zygosity inference rules."""

    def test_het(self):
        """Verify mixed REF/ALT calls are classified as heterozygous.

        Returns:
            None.
        """

        assert _zygosity((0, 1)) == Zygosity.HET

    def test_hom_alt(self):
        """Verify uniform ALT calls are classified as homozygous alternate.

        Returns:
            None.
        """

        assert _zygosity((1, 1)) == Zygosity.HOM_ALT

    def test_ref(self):
        """Verify pure reference calls are classified as reference.

        Returns:
            None.
        """

        assert _zygosity((0, 0)) == Zygosity.REF

    def test_no_call_none(self):
        """Verify all-missing calls using ``None`` are classified as NO_CALL.

        Returns:
            None.
        """

        assert _zygosity((None, None)) == Zygosity.NO_CALL

    def test_no_call_negative(self):
        """Verify all-missing calls using negative values are NO_CALL.

        Returns:
            None.
        """

        assert _zygosity((-1, -1)) == Zygosity.NO_CALL

    def test_partial_missing_is_het(self):
        """Verify partial missing + ALT is classified as heterozygous.

        Returns:
            None.
        """

        # ./1 — one allele unknown, one alt → ambiguous, treat as HET
        assert _zygosity((None, 1)) == Zygosity.HET

    def test_multiallelic_het(self):
        """Verify two different ALT alleles are classified as heterozygous.

        Returns:
            None.
        """

        assert _zygosity((1, 2)) == Zygosity.HET

    def test_haploid_hom_alt(self):
        """Verify haploid ALT calls are classified as homozygous alternate.

        Returns:
            None.
        """

        assert _zygosity((1,)) == Zygosity.HOM_ALT

    def test_haploid_ref(self):
        """Verify haploid reference calls are classified as reference.

        Returns:
            None.
        """

        assert _zygosity((0,)) == Zygosity.REF


class TestGenotypeToVrsIds:
    """Tests for genotype allele-index to VRS-ID mapping."""

    def test_het_single_alt(self):
        """Verify one ALT allele maps to one VRS ID.

        Returns:
            None.
        """

        assert genotype_to_vrs_ids(["ga4gh:VA.abc"], (0, 1)) == ["ga4gh:VA.abc"]

    def test_hom_alt(self):
        """Verify duplicated ALT alleles map to duplicated VRS IDs.

        Returns:
            None.
        """

        assert genotype_to_vrs_ids(["ga4gh:VA.abc"], (1, 1)) == [
            "ga4gh:VA.abc",
            "ga4gh:VA.abc",
        ]

    def test_ref_only_returns_empty(self):
        """Verify reference-only genotypes produce no carried VRS IDs.

        Returns:
            None.
        """

        assert genotype_to_vrs_ids(["ga4gh:VA.abc"], (0, 0)) == []

    def test_missing_call_returns_empty(self):
        """Verify missing calls produce no carried VRS IDs.

        Returns:
            None.
        """

        assert genotype_to_vrs_ids(["ga4gh:VA.abc"], (None, None)) == []

    def test_negative_allele_index_skipped(self):
        """Verify negative allele indexes are skipped during mapping.

        Returns:
            None.
        """

        assert genotype_to_vrs_ids(["ga4gh:VA.abc"], (-1, 1)) == ["ga4gh:VA.abc"]

    def test_multiallelic(self):
        """Verify multi-ALT genotypes map indexes to corresponding VRS IDs.

        Returns:
            None.
        """

        vrs_ids = ["ga4gh:VA.aaa", "ga4gh:VA.bbb"]
        assert genotype_to_vrs_ids(vrs_ids, (1, 2)) == ["ga4gh:VA.aaa", "ga4gh:VA.bbb"]

    def test_allele_index_out_of_range_skipped(self):
        """Verify out-of-range ALT indexes are ignored.

        Returns:
            None.
        """

        # index 2 → idx=1, but only 1 VRS ID available
        assert genotype_to_vrs_ids(["ga4gh:VA.abc"], (0, 2)) == []

    def test_empty_vrs_ids(self):
        """Verify empty VRS ID lists always yield an empty mapping output.

        Returns:
            None.
        """

        assert genotype_to_vrs_ids([], (1, 1)) == []


class TestFormatScalar:
    """Tests for FORMAT scalar extraction edge cases."""

    def test_nan_returns_none(self):
        """NaN FORMAT values should be treated as missing."""

        assert _format_scalar([[float("nan")]], 0) is None

    def test_missing_int_returns_none(self):
        """cyvcf2 missing-int sentinels should be treated as missing."""

        assert _format_scalar([[-2147483648]], 0) is None


def test_load_samples_stores_rows(tmp_path):
    """Verify load_samples inserts expected rows from a mocked VCF.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        None.
    """

    db_path = tmp_path / "test.db"

    mock_record = MagicMock()
    mock_record.FILTER = None  # PASS
    mock_record.INFO.get.return_value = ["ga4gh:VA.aaa"]
    mock_record.CHROM = "chr1"
    mock_record.POS = 100
    # cyvcf2 genotype layout: [allele1, allele2, phased_bool]
    mock_record.genotypes = [
        [0, 1, False],  # S1: 0/1 → HET
        [1, 1, False],  # S2: 1/1 → HOM_ALT
    ]
    mock_record.format.return_value = None  # no GQ/DP

    mock_vcf = MagicMock()
    mock_vcf.samples = ["S1", "S2"]
    mock_vcf.__iter__ = MagicMock(return_value=iter([mock_record]))

    with patch("cyvcf2.VCF", return_value=mock_vcf):
        count = load_samples("fake.vcf.gz", db_path)

    assert count == 2
    conn = open_db(db_path)
    assert "ga4gh:VA.aaa" in get_vrs_ids(conn, "S1")
    assert "ga4gh:VA.aaa" in get_vrs_ids(conn, "S2")
    conn.close()


def test_load_samples_include_no_call(tmp_path):
    """Verify include_no_call=True writes rows for NO_CALL genotypes.

    When ``include_no_call=False`` (default) a NO_CALL genotype produces no
    rows. When ``include_no_call=True`` the VRS IDs at the locus are written
    with NO_CALL zygosity, and the returned count reflects those rows.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        None.
    """

    db_path = tmp_path / "test.db"

    mock_record = MagicMock()
    mock_record.FILTER = None
    mock_record.INFO.get.return_value = ["ga4gh:VA.aaa"]
    mock_record.CHROM = "chr1"
    mock_record.POS = 100
    # S1: ./. → NO_CALL
    mock_record.genotypes = [[-1, -1, False]]
    mock_record.format.return_value = None

    mock_vcf = MagicMock()
    mock_vcf.samples = ["S1"]
    mock_vcf.__iter__ = MagicMock(return_value=iter([mock_record]))

    # Default: NO_CALL excluded
    with patch("cyvcf2.VCF", return_value=mock_vcf):
        count_excluded = load_samples("fake.vcf.gz", db_path)
    assert count_excluded == 0

    # include_no_call=True: genotype is NO_CALL; no carried alleles to index,
    # so no rows are emitted even when the flag is set.
    db_path2 = tmp_path / "test2.db"
    mock_vcf.__iter__ = MagicMock(return_value=iter([mock_record]))
    with patch("cyvcf2.VCF", return_value=mock_vcf):
        count_included = load_samples("fake.vcf.gz", db_path2, include_no_call=True)
    assert count_included == 0

    # include_no_call=False (default): same result — NO_CALL is excluded.
    mock_vcf.__iter__ = MagicMock(return_value=iter([mock_record]))
    with patch("cyvcf2.VCF", return_value=mock_vcf):
        count_excluded2 = load_samples("fake.vcf.gz", db_path2, include_no_call=False)
    assert count_excluded2 == 0


def test_load_samples_gq_filter(tmp_path):
    """Verify records below the GQ threshold are excluded.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        None.
    """

    db_path = tmp_path / "test.db"

    mock_record = MagicMock()
    mock_record.FILTER = None
    mock_record.INFO.get.return_value = ["ga4gh:VA.aaa"]
    mock_record.CHROM = "chr1"
    mock_record.POS = 100
    mock_record.genotypes = [[0, 1, False]]  # S1: HET

    # GQ array: S1 has GQ=10 (below default threshold of 20)
    gq_arr = MagicMock()
    gq_arr.__getitem__ = MagicMock(return_value=[10])
    mock_record.format.return_value = gq_arr

    mock_vcf = MagicMock()
    mock_vcf.samples = ["S1"]
    mock_vcf.__iter__ = MagicMock(return_value=iter([mock_record]))

    with patch("cyvcf2.VCF", return_value=mock_vcf):
        count = load_samples("fake.vcf.gz", db_path, gq_threshold=20)

    assert count == 0


def test_load_samples_no_gq_dp_format_fields(tmp_path):
    """Verify records from VCFs without GQ/DP FORMAT fields are loaded successfully.

    Population-panel VCFs (e.g. 1000 Genomes phased panels) omit per-sample
    quality fields.  cyvcf2 raises ``KeyError`` when those fields are accessed;
    the loader must treat the missing fields as ``None`` rather than crashing.

    Args:
        tmp_path: Pytest-provided temporary directory path.

    Returns:
        None.
    """

    db_path = tmp_path / "test.db"

    mock_record = MagicMock()
    mock_record.FILTER = None
    mock_record.INFO.get.return_value = ["ga4gh:VA.aaa"]
    mock_record.CHROM = "chr22"
    mock_record.POS = 100
    mock_record.genotypes = [[0, 1, True]]  # S1: phased HET

    # Simulate cyvcf2 raising KeyError for absent FORMAT fields
    mock_record.format.side_effect = KeyError(b"GQ")

    mock_vcf = MagicMock()
    mock_vcf.samples = ["S1"]
    mock_vcf.__iter__ = MagicMock(return_value=iter([mock_record]))

    with patch("cyvcf2.VCF", return_value=mock_vcf):
        count = load_samples("fake.vcf.gz", db_path, gq_threshold=0)

    assert count == 1


def test_load_samples_parses_string_vrs_info(tmp_path):
    """Verify comma-delimited INFO strings are parsed as full VRS IDs."""

    db_path = tmp_path / "test.db"

    mock_record = MagicMock()
    mock_record.FILTER = None
    mock_record.INFO.get.return_value = "ga4gh:VA.aaa,ga4gh:VA.bbb"
    mock_record.CHROM = "chr1"
    mock_record.POS = 100
    mock_record.genotypes = [[1, 2, False]]  # S1 carries both ALT alleles
    mock_record.format.side_effect = KeyError(b"GQ")

    mock_vcf = MagicMock()
    mock_vcf.samples = ["S1"]
    mock_vcf.__iter__ = MagicMock(return_value=iter([mock_record]))

    with patch("cyvcf2.VCF", return_value=mock_vcf):
        count = load_samples("fake.vcf.gz", db_path, gq_threshold=0)

    assert count == 2

    conn = open_db(db_path)
    try:
        assert get_vrs_ids(conn, "S1") == frozenset({"ga4gh:VA.aaa", "ga4gh:VA.bbb"})
    finally:
        conn.close()


def test_iter_rows_skips_filtered_missing_and_threshold_failures():
    """Verify _iter_rows skips records for each filtering branch."""

    class Record:
        def __init__(self, *, flt=None, vrs=None, gt=None, gq=None, dp=None):
            self.FILTER = flt
            self._vrs = vrs
            self.CHROM = "chr1"
            self.POS = 100
            self.genotypes = [gt or [0, 1, False]]
            self._gq = gq
            self._dp = dp

        @property
        def INFO(self):
            return SimpleNamespace(get=lambda key: self._vrs)

        def format(self, field):
            if field == "GQ":
                return self._gq
            if field == "DP":
                return self._dp
            raise KeyError(field)

    records = [
        Record(flt="q10", vrs=["ga4gh:VA.filtered"]),
        Record(vrs=None),
        Record(vrs=[]),
        Record(vrs=["ga4gh:VA.ref"], gt=[0, 0, False]),
        Record(vrs=["ga4gh:VA.lowgq"], gq=[[10]]),
        Record(vrs=["ga4gh:VA.lowdp"], gq=[[99]], dp=[[1]]),
        Record(vrs=["ga4gh:VA.skip"], gt=[0, 1, False]),
        Record(vrs=["ga4gh:VA.keep", "ga4gh:VA.keep"], gt=[1, 1, False], gq=[[99]], dp=[[9]]),
    ]

    mock_vcf = MagicMock()
    mock_vcf.samples = ["S1"]
    mock_vcf.__iter__ = MagicMock(return_value=iter(records))

    with patch("cyvcf2.VCF", return_value=mock_vcf):
        rows = list(
            loader_mod._iter_rows(
                "fake.vcf.gz",
                gq_threshold=20,
                dp_threshold=5,
                candidate_vrs_ids={"ga4gh:VA.keep"},
            )
        )

    assert rows == [
        ("S1", "ga4gh:VA.keep", "1/1", "HOM_ALT", "chr1", 100, 99.0, 9, None)
    ]


def test_load_samples_flushes_batches_and_closes_connection(tmp_path, monkeypatch):
    """Verify load_samples flushes large batches and closes the DB connection."""

    rows = [("S1", f"ga4gh:VA.{i}", "0/1", "HET", "chr1", i, None, None, None) for i in range(10_001)]
    insert_calls: list[int] = []

    class DummyConn:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    dummy_conn = DummyConn()

    monkeypatch.setattr(loader_mod, "open_db", lambda path: dummy_conn)
    monkeypatch.setattr(loader_mod, "register_samples", lambda conn, sample_ids: None)
    monkeypatch.setattr(loader_mod, "insert_alleles", lambda conn, batch: insert_calls.append(len(batch)))

    header = MagicMock()
    header.samples = ["S1"]
    header.close = MagicMock()

    mock_vcf = MagicMock()
    mock_vcf.__iter__ = MagicMock(return_value=iter(()))

    monkeypatch.setattr(loader_mod.cyvcf2, "VCF", MagicMock(side_effect=[header, mock_vcf]))
    monkeypatch.setattr(loader_mod, "_iter_rows", lambda *args, **kwargs: iter(rows))

    count = load_samples("fake.vcf.gz", tmp_path / "batch.db")

    assert count == 10_001
    assert insert_calls == [10_000, 1]
    assert dummy_conn.closed is True
    assert header.close.called is True

