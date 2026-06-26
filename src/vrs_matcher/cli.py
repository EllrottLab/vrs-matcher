"""Command-line interface for loading and matching samples.

The commands in this module provide operational access to the ingestion and
matching layers.
"""

import click

from .db import open_db
from .loader import DEFAULT_DP_THRESHOLD, DEFAULT_GQ_THRESHOLD, load_samples
from .matcher import match_against_all, match_pair
from .plugins import PluginError, list_plugins


@click.group()
def cli() -> None:
    """Register the root CLI command group.

    Returns:
        None.
    """


@cli.command("load-samples")
@click.argument("vcf", type=click.Path(exists=True))
@click.option("--db", required=True, type=click.Path(), help="Path to SQLite database.")
@click.option("--source-dataset", default=None, help="Label for the source dataset.")
@click.option(
    "--gq",
    default=DEFAULT_GQ_THRESHOLD,
    show_default=True,
    help="Minimum genotype quality (GQ) to include.",
)
@click.option(
    "--dp",
    default=DEFAULT_DP_THRESHOLD,
    show_default=True,
    help="Minimum read depth (DP) to include.",
)
def load_samples_cmd(vcf: str, db: str, source_dataset: str, gq: float, dp: int) -> None:
    """Load a VRS-annotated VCF into the sample-allele index.

    Args:
        vcf: Path to VCF input file.
        db: Path to SQLite database file.
        source_dataset: Optional source dataset label.
        gq: Minimum genotype quality threshold.
        dp: Minimum depth threshold.

    Returns:
        None.
    """

    n = load_samples(vcf, db, source_dataset=source_dataset, gq_threshold=gq, dp_threshold=dp)
    click.echo(f"Loaded {n} allele records into {db}")


@cli.command("match-samples")
@click.argument("sample_a")
@click.argument("sample_b")
@click.option("--db", required=True, type=click.Path(exists=True), help="Path to SQLite database.")
@click.option(
    "--algorithm",
    default="identity",
    show_default=True,
    help="Matching plugin algorithm name.",
)
@click.option(
    "--plugin-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to a local Python plugin file with create_plugin().",
)
def match_samples_cmd(
    sample_a: str,
    sample_b: str,
    db: str,
    algorithm: str,
    plugin_file: str | None,
) -> None:
    """Compute and print similarity metrics for two samples.

    Args:
        sample_a: First sample identifier.
        sample_b: Second sample identifier.
        db: Path to SQLite database file.

    Returns:
        None.
    """

    conn = open_db(db)
    try:
        result = match_pair(conn, sample_a, sample_b, algorithm=algorithm, plugin_file=plugin_file)
    except PluginError as exc:
        raise click.ClickException(str(exc)) from exc
    except KeyError as exc:
        raise click.ClickException(f"Sample not found in index: {exc.args[0]}") from exc
    finally:
        conn.close()

    click.echo(f"Jaccard:              {result.jaccard:.4f}")
    click.echo(f"Weighted concordance: {result.weighted_concordance:.4f}")
    click.echo(f"Shared variants:      {len(result.shared_vrs_ids)}")
    click.echo(f"Total alleles (A/B):  {result.total_a} / {result.total_b}")


@cli.command("match-sample")
@click.argument("sample_id")
@click.option(
    "--against",
    default="all",
    show_default=True,
    type=click.Choice(["all"]),
    help="Comparison target (currently only 'all' is supported).",
)
@click.option(
    "--top",
    default=20,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of top matches to return.",
)
@click.option("--db", required=True, type=click.Path(exists=True), help="Path to SQLite database.")
@click.option(
    "--algorithm",
    default="identity",
    show_default=True,
    help="Matching plugin algorithm name.",
)
@click.option(
    "--plugin-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to a local Python plugin file with create_plugin().",
)
def match_sample_cmd(
    sample_id: str,
    against: str,
    top: int,
    db: str,
    algorithm: str,
    plugin_file: str | None,
) -> None:
    """Match one sample against all others and print ranked results.

    Args:
        sample_id: Query sample identifier.
        against: Target scope selector (currently only ``"all"``).
        top: Maximum number of matches to display.
        db: Path to SQLite database file.

    Returns:
        None.
    """

    conn = open_db(db)
    try:
        results = match_against_all(
            conn,
            sample_id,
            top_n=top,
            algorithm=algorithm,
            plugin_file=plugin_file,
        )
    except PluginError as exc:
        raise click.ClickException(str(exc)) from exc
    except KeyError as exc:
        raise click.ClickException(f"Sample not found in index: {exc.args[0]}") from exc
    finally:
        conn.close()

    if not results:
        click.echo("No other samples found in index.")
        return

    click.echo(f"{'Sample':<30} {'Jaccard':>8} {'WConc':>8} {'Shared':>8}")
    click.echo("-" * 60)
    for r in results:
        other = r.sample_b if r.sample_a == sample_id else r.sample_a
        click.echo(
            f"{other:<30} {r.jaccard:>8.4f} {r.weighted_concordance:>8.4f}"
            f" {len(r.shared_vrs_ids):>8}"
        )


@cli.command("shared-variants")
@click.argument("sample_a")
@click.argument("sample_b")
@click.option("--db", required=True, type=click.Path(exists=True), help="Path to SQLite database.")
@click.option(
    "--algorithm",
    default="identity",
    show_default=True,
    help="Matching plugin algorithm name.",
)
@click.option(
    "--plugin-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Path to a local Python plugin file with create_plugin().",
)
def shared_variants_cmd(
    sample_a: str,
    sample_b: str,
    db: str,
    algorithm: str,
    plugin_file: str | None,
) -> None:
    """Print shared VRS IDs for two samples.

    Args:
        sample_a: First sample identifier.
        sample_b: Second sample identifier.
        db: Path to SQLite database file.

    Returns:
        None.
    """

    conn = open_db(db)
    try:
        result = match_pair(conn, sample_a, sample_b, algorithm=algorithm, plugin_file=plugin_file)
    except PluginError as exc:
        raise click.ClickException(str(exc)) from exc
    except KeyError as exc:
        raise click.ClickException(f"Sample not found in index: {exc.args[0]}") from exc
    finally:
        conn.close()

    for vrs_id in sorted(result.shared_vrs_ids):
        click.echo(vrs_id)


@cli.group("plugins")
def plugins_cmd() -> None:
    """Commands for inspecting discovered matcher plugins."""


@plugins_cmd.command("list")
def plugins_list_cmd() -> None:
    """List available matcher plugin names."""

    for name in list_plugins():
        click.echo(name)
