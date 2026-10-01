"""Explicit phenotype commands; imported only by opt-in CLI paths."""

import json
import sqlite3
from contextlib import closing

import click

from .phenotypes import Phenotypes, export_pairs, import_snapshot, snapshot_diff
from .plugins import PluginError


def emit(value):
    click.echo(json.dumps(value, default=sorted, allow_nan=False, indent=2))


def checked(function, *args, **kwargs):
    try:
        return function(*args, **kwargs)
    except (ValueError, KeyError, OSError, sqlite3.Error, PluginError) as exc:
        raise click.ClickException(str(exc)) from exc


def register(cli):
    @cli.command("stage-gregor-inputs")
    @click.option("--manifest", required=True, type=click.Path(exists=True, dir_okay=False))
    @click.option("--output", required=True, type=click.Path())
    def stage_cmd(manifest, output):
        """Verify and atomically stage already downloaded files (no network)."""
        from .gregor_inputs import stage_inputs

        emit(checked(stage_inputs, manifest, output))

    @cli.command("import-phenotypes")
    @click.option("--db", required=True, type=click.Path(exists=True, dir_okay=False))
    @click.option("--phenotype-db", required=True, type=click.Path(dir_okay=False))
    @click.option("--manifest", required=True, type=click.Path(exists=True, dir_okay=False))
    def import_cmd(db, phenotype_db, manifest):
        """Validate and atomically publish an offline clinical bundle."""
        emit(checked(import_snapshot, db, phenotype_db, manifest))

    @cli.command("phenotype-diff")
    @click.option("--phenotype-db", required=True, type=click.Path(exists=True, dir_okay=False))
    @click.argument("before")
    @click.argument("after")
    def diff_cmd(phenotype_db, before, after):
        """Compare clinical content and mappings across snapshots."""
        emit(checked(snapshot_diff, phenotype_db, before, after))

    @cli.command("phenotype-cohort")
    @click.option("--phenotype-db", required=True, type=click.Path(exists=True, dir_okay=False))
    @click.option("--phenotype-snapshot", required=True)
    @click.option("--ontology", required=True)
    @click.option("--term-id", required=True)
    @click.option(
        "--presence", default="Present", type=click.Choice(["Present", "Absent", "Unknown"])
    )
    @click.option("--exclude-conflicts", is_flag=True)
    def cohort_cmd(phenotype_db, phenotype_snapshot, **kwargs):
        """Select an exact-term clinical cohort without genetic scoring."""
        with closing(checked(Phenotypes, phenotype_db, phenotype_snapshot)) as reader:
            emit(checked(reader.cohort, **kwargs))

    @cli.command("export-phenotype-pairs")
    @click.option("--db", required=True, type=click.Path(exists=True, dir_okay=False))
    @click.option("--phenotype-db", required=True, type=click.Path(exists=True, dir_okay=False))
    @click.option("--phenotype-snapshot", required=True)
    @click.option("--output", required=True, type=click.Path())
    @click.option(
        "--algorithm", default="king-robust", type=click.Choice(["identity", "king-robust"])
    )
    def export_cmd(db, phenotype_db, phenotype_snapshot, output, algorithm):
        """Stream pairs and participant annotations into a new report directory."""
        emit(
            checked(export_pairs, db, phenotype_db, phenotype_snapshot, output, algorithm=algorithm)
        )
