"""Explicitly enabled local GREGoR validation; never fetches protected inputs."""

import json
import os
from contextlib import closing
from pathlib import Path

import pytest

from vrs_matcher.phenotypes import Phenotypes, import_snapshot


@pytest.mark.integration
def test_authorized_gregor_bundle(tmp_path):
    if os.environ.get("VRS_MATCHER_GREGOR_INTEGRATION") != "1":
        pytest.skip("Set VRS_MATCHER_GREGOR_INTEGRATION=1 for authorized local validation")
    # Once opted in, missing configuration/data is a failure, not a successful skip.
    manifest = Path(os.environ["VRS_MATCHER_GREGOR_MANIFEST"]).resolve()
    db = Path(os.environ["VRS_MATCHER_GREGOR_DB"]).resolve()
    truth_path = Path(os.environ["VRS_MATCHER_GREGOR_EXPECTED"]).resolve()
    truth = json.loads(truth_path.read_text())
    sidecar = tmp_path / "clinical.db"
    result = import_snapshot(db, sidecar, manifest)
    with closing(Phenotypes(sidecar, result["snapshot_id"], db)) as reader:
        report = reader.enrich(truth["sample_ids"])
        assert report["observations"] == truth["observations"]
        assert report["counts"] == truth["counts"]
        assert reader.cohort(**truth["predicate"])["sample_ids"] == truth["cohort_sample_ids"]
