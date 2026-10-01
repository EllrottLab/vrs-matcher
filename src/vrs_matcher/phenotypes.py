"""Offline GREGoR snapshots. Clinical reads never enter genetic scoring loops."""

import csv
import hashlib
import json
import math
import sqlite3
from collections import Counter
from contextlib import closing
from datetime import UTC, datetime
from itertools import batched
from pathlib import Path

from .genotypes import file_sha256
from .plugins import PluginContext

_SCHEMA = """
CREATE TABLE IF NOT EXISTS phenotype_snapshot (
 snapshot_id TEXT PRIMARY KEY NOT NULL, manifest TEXT NOT NULL,
 imported_at TEXT NOT NULL, importer_version INTEGER NOT NULL CHECK(importer_version=1),
 binding_stamp TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS phenotype_participant (
 snapshot_id TEXT NOT NULL REFERENCES phenotype_snapshot(snapshot_id),
 source_id TEXT NOT NULL, participant_id TEXT NOT NULL, person_key TEXT NOT NULL,
 family_id TEXT, fields TEXT NOT NULL,
 PRIMARY KEY(snapshot_id,source_id,participant_id)
);
CREATE INDEX IF NOT EXISTS participant_person
 ON phenotype_participant(snapshot_id,person_key);
CREATE TABLE IF NOT EXISTS observation_participant (
 snapshot_id TEXT NOT NULL, sample_id TEXT NOT NULL,
 source_id TEXT NOT NULL, participant_id TEXT NOT NULL,
 PRIMARY KEY(snapshot_id,sample_id),
 FOREIGN KEY(snapshot_id,source_id,participant_id)
 REFERENCES phenotype_participant(snapshot_id,source_id,participant_id)
);
CREATE INDEX IF NOT EXISTS observation_person
 ON observation_participant(snapshot_id,source_id,participant_id,sample_id);
CREATE TABLE IF NOT EXISTS phenotype_assertion (
 snapshot_id TEXT NOT NULL, assertion_id TEXT NOT NULL, source_id TEXT NOT NULL,
 participant_id TEXT NOT NULL, ontology TEXT NOT NULL, term_id TEXT NOT NULL,
 presence TEXT NOT NULL CHECK(presence IN ('Present','Absent','Unknown')),
 phenotype_id TEXT, fields TEXT NOT NULL, file_sha256 TEXT NOT NULL, row_number INTEGER NOT NULL,
 PRIMARY KEY(snapshot_id,assertion_id),
 FOREIGN KEY(snapshot_id,source_id,participant_id)
 REFERENCES phenotype_participant(snapshot_id,source_id,participant_id)
);
CREATE INDEX IF NOT EXISTS assertion_participant
 ON phenotype_assertion(snapshot_id,source_id,participant_id);
CREATE INDEX IF NOT EXISTS assertion_term
 ON phenotype_assertion(snapshot_id,ontology,term_id,presence,source_id,participant_id);
"""
_PREFIX = {
    "HPO": "HP:",
    "MONDO": "MONDO:",
    "OMIM": "OMIM:",
    "ORPHANET": "ORPHA:",
    "SNOMED": "SNOMED:",
    "ICD10": "ICD10:",
}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def read_only(path):
    conn = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def stamp(path):
    st = Path(path).stat()
    return [st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns]


def genetic_binding(path, artifact_id):
    """Hash a quiescent, checkpointed index once, never a live WAL database."""
    if not artifact_id or Path(str(path) + "-wal").exists():
        raise ValueError("A named, checkpointed genetic artifact is required.")
    before = stamp(path)
    sha = file_sha256(path)
    with closing(read_only(path)) as conn:
        registry = digest([r[0] for r in conn.execute("SELECT sample_id FROM samples ORDER BY 1")])
    if before != stamp(path):
        raise ValueError("Genetic artifact changed during verification.")
    return {"artifact_id": artifact_id, "sha256": sha, "registry_sha256": registry}, before


def contract(version):
    if version not in ("1.11", "1.12"):
        raise ValueError("Supported GREGoR model versions: 1.11, 1.12.")
    path = Path(__file__).parent / "schemas" / f"gregor-{version}.json"
    return json.loads(path.read_text())


def rows(path, required, columns=None):
    """Read deposition or scalar-reference Terra TSVs without guessing identifiers."""
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        headers = reader.fieldnames or []
        names = [h.removeprefix("entity:") for h in headers]
        if len(set(names)) != len(names) or not required <= set(names):
            raise ValueError(f"{Path(path).name}: missing/duplicate columns")
        for number, raw in enumerate(reader, 2):
            if None in raw or any(v is None for v in raw.values()):
                raise ValueError(f"{Path(path).name}:{number}: malformed TSV row")
            row = dict(zip(names, raw.values(), strict=True))
            for key in ("participant_id",):
                value = row.get(key, "")
                if value.startswith("{"):
                    try:
                        ref = json.loads(value)
                        if (
                            set(ref) != {"entityType", "entityName"}
                            or ref["entityType"] != "participant"
                        ):
                            raise ValueError
                        row[key] = ref["entityName"]
                    except (ValueError, TypeError):
                        raise ValueError(
                            f"{Path(path).name}:{number}: invalid entity reference"
                        ) from None
            if any(
                not isinstance(row[k], str) or not row[k] or row[k] != row[k].strip()
                for k in required
            ):
                raise ValueError(f"{Path(path).name}:{number}: missing/invalid required value")
            if columns:
                for col in columns:
                    key, value = col["column"], row.get(col["column"], "")
                    if not value:
                        continue
                    values = value.split(col.get("multi_value_delimiter", "\0"))
                    if "enumerations" in col and any(v not in col["enumerations"] for v in values):
                        raise ValueError(f"{Path(path).name}:{number}: invalid {key}")
                    if col.get("data_type") == "float":
                        try:
                            if not math.isfinite(float(value)) or float(value) < 0:
                                raise ValueError
                        except ValueError:
                            raise ValueError(f"{Path(path).name}:{number}: invalid {key}") from None
            yield number, row


def import_snapshot(db, phenotype_db, manifest_path):
    """Import an entire clinical bundle atomically; the genetic artifact is read-only."""
    if Path(db).resolve() == Path(phenotype_db).resolve():
        raise ValueError("Phenotype sidecar must differ from genetic database.")
    root = Path(manifest_path).resolve().parent
    manifest = json.loads(Path(manifest_path).read_text())
    schema = contract(manifest["model_version"])
    binding, verified_stamp = genetic_binding(db, manifest["genetic"]["artifact_id"])
    if manifest["genetic"].get("sha256") != binding["sha256"]:
        raise ValueError("Genetic artifact checksum mismatch.")
    sources = manifest["sources"]
    if not sources or len({s["source_id"] for s in sources}) != len(sources):
        raise ValueError("Source namespaces must be nonempty and unique.")
    paths, published = {}, []
    for source in sources:
        if not isinstance(source["source_id"], str) or not source["source_id"].strip():
            raise ValueError("Empty source namespace.")
        entry = {"source_id": source["source_id"]}
        for table in ("participant", "phenotype", "mapping", "reconciliation"):
            if table == "reconciliation" and table not in source:
                continue
            spec = source[table]
            path = root / spec["path"]
            actual = file_sha256(path)
            if actual != spec["sha256"]:
                raise ValueError(f"{source['source_id']}/{table}: input checksum mismatch")
            paths[source["source_id"], table] = (path, actual)
            entry[table] = {"sha256": actual}
        published.append(entry)
    metadata = {
        "model_version": schema["version"],
        "schema_sha256": schema["upstream_sha256"],
        "canonical_version": 1,
        "genetic": binding,
        "sources": sorted(published, key=lambda s: s["source_id"]),
        "provenance": manifest.get("provenance", {}),
    }
    snapshot = digest(metadata)
    with closing(read_only(db)) as genetic, closing(sqlite3.connect(phenotype_db)) as conn:
        conn.execute("PRAGMA foreign_keys=ON")
        conn.executescript(_SCHEMA)
        conn.execute("BEGIN IMMEDIATE")
        with conn:
            existing = conn.execute(
                "SELECT 1 FROM phenotype_snapshot WHERE snapshot_id=?", (snapshot,)
            ).fetchone()
            if not existing:
                conn.execute(
                    "INSERT INTO phenotype_snapshot VALUES (?,?,?,?,?)",
                    (
                        snapshot,
                        canonical(metadata),
                        datetime.now(UTC).isoformat(),
                        1,
                        canonical(verified_stamp),
                    ),
                )
                for source in published:
                    sid = source["source_id"]
                    cols = schema["tables"]["participant"]
                    # Participant context is a projection, not a full GREGoR deposition validator.
                    for _, row in rows(paths[sid, "participant"][0], {"participant_id"}, cols):
                        pid = row["participant_id"]
                        conn.execute(
                            "INSERT INTO phenotype_participant VALUES (?,?,?,?,?,?)",
                            (
                                snapshot,
                                sid,
                                pid,
                                canonical(["source", sid, pid]),
                                row.get("family_id") or None,
                                canonical(row),
                            ),
                        )
                    if (sid, "reconciliation") in paths:
                        seen = set()
                        for _, row in rows(
                            paths[sid, "reconciliation"][0], {"participant_id", "person_key"}
                        ):
                            if row["participant_id"] in seen:
                                raise ValueError("Duplicate reconciliation assignment.")
                            seen.add(row["participant_id"])
                            result = conn.execute(
                                "UPDATE phenotype_participant SET person_key=? "
                                "WHERE snapshot_id=? AND source_id=? "
                                "AND participant_id=?",
                                (
                                    canonical(["reconciled", row["person_key"]]),
                                    snapshot,
                                    sid,
                                    row["participant_id"],
                                ),
                            )
                            if result.rowcount != 1:
                                raise ValueError("Unknown reconciliation participant.")
                    cols = schema["tables"]["phenotype"]
                    required = {c["column"] for c in cols if c.get("required")}
                    path, sha = paths[sid, "phenotype"]
                    for number, row in rows(path, required, cols):
                        if (
                            not row["term_id"].startswith(_PREFIX[row["ontology"]])
                            or not row["term_id"].split(":", 1)[1]
                        ):
                            raise ValueError(
                                f"{sid}/phenotype:{number}: ontology namespace mismatch"
                            )
                        conn.execute(
                            "INSERT INTO phenotype_assertion VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                            (
                                snapshot,
                                digest([sid, sha, number]),
                                sid,
                                row["participant_id"],
                                row["ontology"],
                                row["term_id"],
                                row["presence"],
                                row.get("phenotype_id"),
                                canonical(row),
                                sha,
                                number,
                            ),
                        )
                    for _, row in rows(paths[sid, "mapping"][0], {"sample_id", "participant_id"}):
                        if not genetic.execute(
                            "SELECT 1 FROM samples WHERE sample_id=?", (row["sample_id"],)
                        ).fetchone():
                            raise ValueError("Mapping references an unregistered observation.")
                        conn.execute(
                            "INSERT INTO observation_participant VALUES (?,?,?,?)",
                            (snapshot, row["sample_id"], sid, row["participant_id"]),
                        )
            for path, sha in paths.values():
                if file_sha256(path) != sha:
                    raise ValueError("Input changed during import.")
            if stamp(db) != verified_stamp:
                raise ValueError("Genetic artifact changed during import.")
        total = genetic.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        counts = {
            t: conn.execute(
                f"SELECT COUNT(*) FROM {t} WHERE snapshot_id=?", (snapshot,)
            ).fetchone()[0]
            for t in ("phenotype_participant", "phenotype_assertion", "observation_participant")
        }
        counts["unmapped_observations"] = total - counts["observation_participant"]
        counts["presence"] = dict(
            conn.execute(
                "SELECT presence, COUNT(*) FROM phenotype_assertion "
                "WHERE snapshot_id=? GROUP BY presence",
                (snapshot,),
            )
        )
        counts["conflicting_terms"] = conn.execute(
            """SELECT COUNT(*) FROM (
          SELECT p.person_key,a.ontology,a.term_id FROM phenotype_assertion a
          JOIN phenotype_participant p USING(snapshot_id,source_id,participant_id)
          WHERE a.snapshot_id=? GROUP BY p.person_key,a.ontology,a.term_id
          HAVING COUNT(DISTINCT a.presence)>1)""",
            (snapshot,),
        ).fetchone()[0]
        counts["repeated_assertions"] = conn.execute(
            """SELECT COALESCE(SUM(n-1),0) FROM (
          SELECT COUNT(*) n FROM phenotype_assertion WHERE snapshot_id=?
          GROUP BY source_id,participant_id,json_remove(fields,'$.phenotype_id'))""",
            (snapshot,),
        ).fetchone()[0]
        return {"snapshot_id": snapshot, "already_present": bool(existing), **counts}


class Phenotypes:
    """Explicit snapshot reader; owns a read-only sidecar connection."""

    def __init__(self, path, snapshot, db=None):
        self.conn = read_only(path)
        self.snapshot = snapshot
        try:
            row = self.conn.execute(
                "SELECT * FROM phenotype_snapshot WHERE snapshot_id=?", (snapshot,)
            ).fetchone()
            if row is None:
                raise ValueError("Unknown phenotype snapshot.")
            self.manifest = json.loads(row["manifest"])
            if db is not None:
                # Unchanged local immutable artifact: O(1). Relocation/change: verify bytes.
                if Path(str(db) + "-wal").exists():
                    raise ValueError("Use a checkpointed genetic artifact.")
                if stamp(db) != json.loads(row["binding_stamp"]):
                    binding, _ = genetic_binding(db, self.manifest["genetic"]["artifact_id"])
                    if binding != self.manifest["genetic"]:
                        raise ValueError("Genetic artifact does not match phenotype snapshot.")
        except BaseException:
            self.conn.close()
            raise

    def close(self):
        self.conn.close()

    def _batches(self, ids):
        size = min(500, self.conn.getlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER) - 1)
        if size < 1:
            raise ValueError("SQLite parameter limit too small.")
        return batched(ids, size, strict=False)

    def people(self, keys):
        result = {}
        for batch in self._batches(keys):
            marks = ",".join("?" for _ in batch)
            for row in self.conn.execute(
                "SELECT * FROM phenotype_participant "
                f"WHERE snapshot_id=? AND person_key IN ({marks}) "
                "ORDER BY source_id,participant_id",
                (self.snapshot, *batch),
            ):
                person = result.setdefault(row["person_key"], {"records": [], "assertions": []})
                person["records"].append(
                    {
                        "source_id": row["source_id"],
                        "participant_id": row["participant_id"],
                        "fields": json.loads(row["fields"]),
                    }
                )
            for row in self.conn.execute(
                """SELECT a.*,p.person_key FROM phenotype_participant p
              JOIN phenotype_assertion a USING(snapshot_id,source_id,participant_id)
              WHERE p.snapshot_id=? AND p.person_key IN ("""
                + marks
                + """)
              ORDER BY p.person_key,a.source_id,a.participant_id,
                       a.ontology,a.term_id,a.assertion_id""",
                (self.snapshot, *batch),
            ):
                data = dict(row)
                key = data.pop("person_key")
                data["fields"] = json.loads(data["fields"])
                result[key]["assertions"].append(data)
        for person in result.values():
            statuses = {}
            for a in person["assertions"]:
                statuses.setdefault((a["ontology"], a["term_id"]), set()).add(a["presence"])
            person["terms"] = [
                {"ontology": o, "term_id": t, "statuses": sorted(v), "conflict": len(v) > 1}
                for (o, t), v in sorted(statuses.items())
            ]
            person["status"] = "recorded" if statuses else "no_assertions"
        return result

    def enrich(self, sample_ids):
        mappings = dict.fromkeys(sorted(set(sample_ids)))
        for batch in self._batches(mappings):
            marks = ",".join("?" for _ in batch)
            for row in self.conn.execute(
                """SELECT m.sample_id,p.person_key
              FROM observation_participant m JOIN phenotype_participant p
              USING(snapshot_id,source_id,participant_id)
              WHERE m.snapshot_id=? AND m.sample_id IN ("""
                + marks
                + ")",
                (self.snapshot, *batch),
            ):
                mappings[row["sample_id"]] = row["person_key"]
        people = self.people(sorted({p for p in mappings.values() if p is not None}))
        families = {
            (r["source_id"], r["fields"]["family_id"])
            for p in people.values()
            for r in p["records"]
            if r["fields"].get("family_id")
        }
        return {
            "snapshot_id": self.snapshot,
            "manifest": self.manifest,
            "observations": {
                s: {"person_key": p, "status": "mapped" if p else "unmapped"}
                for s, p in mappings.items()
            },
            "participants": people,
            "counts": {
                "people": len(people),
                "scoped_families": len(families),
                "people_without_family": sum(
                    not any(r["fields"].get("family_id") for r in p["records"])
                    for p in people.values()
                ),
                "unmapped_observations": sum(p is None for p in mappings.values()),
            },
        }

    def cohort(self, ontology, term_id, presence="Present", exclude_conflicts=False):
        if ontology not in _PREFIX or presence not in ("Present", "Absent", "Unknown"):
            raise ValueError("Invalid cohort ontology or presence.")
        statuses = {}
        for row in self.conn.execute(
            """SELECT DISTINCT p.person_key,a.presence
          FROM phenotype_assertion a JOIN phenotype_participant p
          USING(snapshot_id,source_id,participant_id)
          WHERE a.snapshot_id=? AND a.ontology=? AND a.term_id=?""",
            (self.snapshot, ontology, term_id),
        ):
            statuses.setdefault(row[0], set()).add(row[1])
        selected = {
            p
            for p, values in statuses.items()
            if presence in values and (not exclude_conflicts or len(values) == 1)
        }
        ids = []
        for batch in self._batches(sorted(selected)):
            marks = ",".join("?" for _ in batch)
            ids.extend(
                r[0]
                for r in self.conn.execute(
                    """SELECT m.sample_id
              FROM phenotype_participant p JOIN observation_participant m
              USING(snapshot_id,source_id,participant_id)
              WHERE p.snapshot_id=? AND p.person_key IN ("""
                    + marks
                    + ")",
                    (self.snapshot, *batch),
                )
            )
        counts = Counter(
            next(iter(v)) + "_only" if len(v) == 1 else "conflicting" for v in statuses.values()
        )
        total = self.conn.execute(
            "SELECT COUNT(DISTINCT person_key) FROM phenotype_participant WHERE snapshot_id=?",
            (self.snapshot,),
        ).fetchone()[0]
        counts["unrecorded"] = total - len(statuses)
        ids = sorted(set(ids))
        return {
            "sample_ids": ids,
            "sample_ids_sha256": digest(ids),
            "count": len(ids),
            "selected_people": len(selected),
            "groups": dict(counts),
            "predicate": {
                "ontology": ontology,
                "term_id": term_id,
                "presence": presence,
                "exclude_conflicts": exclude_conflicts,
            },
            "snapshot_id": self.snapshot,
        }


def snapshot_diff(path, before, after):
    """SQL set difference ignores row order, source IDs and duplicate assertions."""
    with closing(Phenotypes(path, before)) as old, closing(Phenotypes(path, after)) as new:
        # Clinical field projection excludes generated source IDs; blank optional values
        # are kept, and mappings/context changes are reported independently.
        query = """SELECT source_id,participant_id,json_remove(fields,'$.phenotype_id')
                   FROM phenotype_assertion WHERE snapshot_id=?"""
        changes = {}
        for name, sql in [
            ("assertions", query),
            (
                "participants",
                "SELECT source_id,participant_id,person_key,fields "
                "FROM phenotype_participant WHERE snapshot_id=?",
            ),
            (
                "mappings",
                "SELECT sample_id,source_id,participant_id "
                "FROM observation_participant WHERE snapshot_id=?",
            ),
        ]:
            changes[name] = {}
            for label, a, b in [("removed", before, after), ("added", after, before)]:
                changes[name][label] = [
                    list(r) for r in old.conn.execute(sql + " EXCEPT " + sql, (a, b))
                ]
        return {"before": old.snapshot, "after": new.snapshot, "changes": changes}


class CohortContext(PluginContext):
    def __init__(self, conn, samples):
        super().__init__(conn)
        self.samples = sorted(set(samples))

    def list_samples(self):
        return self.samples


def match_report(
    db,
    phenotype_db,
    snapshot,
    sample_a,
    sample_b=None,
    *,
    algorithm="identity",
    plugin_file=None,
    top_n=20,
    predicate=None,
    candidate_vrs_ids=None,
):
    """Opt-in orchestration: genetic scoring first, batched clinical enrichment last."""
    from dataclasses import asdict

    from .models import KinshipMatches
    from .plugins import resolve_plugin

    with (
        closing(Phenotypes(phenotype_db, snapshot, db)) as clinical,
        closing(read_only(db)) as conn,
    ):
        if predicate and (
            sample_b is not None or plugin_file or algorithm not in ("identity", "king-robust")
        ):
            raise ValueError("Cohort selection supports built-in one-versus-all matching only.")
        plugin = resolve_plugin(name=algorithm, plugin_file=plugin_file)
        context = PluginContext(conn)
        cohort = None
        if predicate:
            cohort = clinical.cohort(**predicate)
            registered = set(context.list_samples())
            peers = sorted(registered.intersection(cohort["sample_ids"]))
            cohort.update(
                sample_ids=peers,
                count=len(peers),
                sample_ids_sha256=digest(peers),
                query_outside_cohort=sample_a not in peers,
            )
            context = CohortContext(conn, [*peers, sample_a] if sample_a in registered else peers)
        if sample_b is not None:
            result = plugin.match_pair(
                context, sample_a, sample_b, candidate_vrs_ids=candidate_vrs_ids
            )
            payload, ids = asdict(result), [sample_a, sample_b]
        else:
            result = plugin.match_against_all(
                context, sample_a, top_n=top_n, candidate_vrs_ids=candidate_vrs_ids
            )
            all_results = (
                result.matches + result.unscorable if isinstance(result, KinshipMatches) else result
            )
            ids = [sample_a, *(r.sample_b for r in all_results)]
            payload = (
                asdict(result)
                if isinstance(result, KinshipMatches)
                else [asdict(r) for r in result]
            )
        return {
            "report_version": 1,
            "genetic": payload,
            "phenotypes": clinical.enrich(ids),
            "cohort": cohort,
            "parameters": {
                "algorithm": plugin.name,
                "sample_a": sample_a,
                "sample_b": sample_b,
                "top_n": top_n,
                "candidate_vrs_ids": sorted(candidate_vrs_ids)
                if candidate_vrs_ids is not None
                else None,
            },
        }


def export_pairs(db, phenotype_db, snapshot, output, *, algorithm="king-robust", predicate=None):
    """Publish a directory only after both bounded streams and their manifest succeed."""
    import itertools
    import tempfile
    from dataclasses import asdict

    from .plugins import resolve_plugin

    destination = Path(output)
    if destination.exists():
        raise ValueError("Export destination must not already exist.")
    if algorithm not in ("identity", "king-robust"):
        raise ValueError("All-pairs export supports built-in algorithms only.")
    with (
        closing(Phenotypes(phenotype_db, snapshot, db)) as clinical,
        closing(read_only(db)) as conn,
    ):
        context = PluginContext(conn)
        cohort = clinical.cohort(**predicate) if predicate else None
        if cohort:
            context = CohortContext(conn, set(context.list_samples()) & set(cohort["sample_ids"]))
        samples = context.list_samples()
        plugin = resolve_plugin(name=algorithm)
        pairs = (
            plugin.iter_all_pairs(context)
            if algorithm == "king-robust"
            else (plugin.match_pair(context, a, b) for a, b in itertools.combinations(samples, 2))
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=destination.parent) as staging:
            stage = Path(staging) / "report"
            stage.mkdir()
            count = 0
            with (stage / "pairs.jsonl").open("w") as stream:
                for result in pairs:
                    stream.write(json.dumps(asdict(result), default=sorted, allow_nan=False) + "\n")
                    count += 1
            # Use a temporary SQLite set, not an unbounded Python assertion cache.
            clinical.conn.execute("CREATE TEMP TABLE selected_people(person_key TEXT PRIMARY KEY)")
            with (stage / "observations.jsonl").open("w") as stream:
                for batch in clinical._batches(samples):
                    marks = ",".join("?" for _ in batch)
                    mappings = dict.fromkeys(batch)
                    for row in clinical.conn.execute(
                        """SELECT m.sample_id,p.person_key
                        FROM observation_participant m JOIN phenotype_participant p
                        USING(snapshot_id,source_id,participant_id)
                        WHERE m.snapshot_id=? AND m.sample_id IN ("""
                        + marks
                        + ")",
                        (snapshot, *batch),
                    ):
                        mappings[row[0]] = row[1]
                    for sample, person in mappings.items():
                        stream.write(canonical({"sample_id": sample, "person_key": person}) + "\n")
                        if person is not None:
                            clinical.conn.execute(
                                "INSERT OR IGNORE INTO selected_people VALUES (?)", (person,)
                            )
            with (stage / "participants.jsonl").open("w") as stream:
                keys = (
                    r[0]
                    for r in clinical.conn.execute(
                        "SELECT person_key FROM selected_people ORDER BY person_key"
                    )
                )
                for batch in clinical._batches(keys):
                    for key, person in clinical.people(batch).items():
                        stream.write(canonical({"person_key": key, **person}) + "\n")
            metadata = {
                "report_version": 1,
                "snapshot_id": snapshot,
                "phenotype_manifest": clinical.manifest,
                "algorithm": algorithm,
                "cohort": cohort,
                "samples": len(samples),
                "pairs": count,
                "files": {p.name: file_sha256(p) for p in sorted(stage.iterdir())},
            }
            (stage / "manifest.json").write_text(canonical(metadata) + "\n")
            stage.rename(destination)
    return metadata
