"""Unit tests for core model types.

These tests validate enum values and dataclass field behavior for model objects
used throughout indexing and matching logic.
"""

from vrs_matcher.models import GenotypeState, SampleGenotype, Zygosity


def test_zygosity_values():
    """Verify expected literal values for all zygosity enum members.

    Returns:
        None.
    """

    assert Zygosity.HET.value == "HET"
    assert Zygosity.HOM_ALT.value == "HOM_ALT"
    assert Zygosity.REF.value == "REF"
    assert Zygosity.NO_CALL.value == "NO_CALL"


def test_zygosity_str_enum_equality():
    """Verify string-based enum construction resolves to the same member.

    Returns:
        None.
    """

    assert Zygosity("HET") == Zygosity.HET


def test_genotype_state_defaults():
    """Verify optional genotype state fields default to ``None``.

    Returns:
        None.
    """

    gs = GenotypeState(gt="0/1", zygosity=Zygosity.HET)
    assert gs.depth is None
    assert gs.gq is None


def test_genotype_state_full():
    """Verify explicit genotype state fields are preserved.

    Returns:
        None.
    """

    gs = GenotypeState(gt="1/1", zygosity=Zygosity.HOM_ALT, gq=40.0, depth=30)
    assert gs.gq == 40.0
    assert gs.depth == 30


def test_sample_genotype_empty():
    """Verify sample genotype starts with an empty allele mapping by default.

    Returns:
        None.
    """

    sg = SampleGenotype(sample_id="S1")
    assert sg.alleles == {}


def test_sample_genotype_with_alleles():
    """Verify sample genotype stores and exposes provided allele states.

    Returns:
        None.
    """

    gs = GenotypeState(gt="1/1", zygosity=Zygosity.HOM_ALT)
    sg = SampleGenotype(sample_id="S1", alleles={"ga4gh:VA.abc": gs})
    assert sg.alleles["ga4gh:VA.abc"].gt == "1/1"
