"""CLI for KEGG module completeness analysis."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .cog_store import get_cog_catalog
from .kegg_api import KeggAPI
from .ko_input import resolve_kos
from .module_completeness import (
    KeggModuleAnalyzer,
    default_tsv_output_path,
    write_module_completeness_tsv,
)

MANIFEST_LARGE_ROW_THRESHOLD = 1000
GENOME_ID_COLUMNS = ("genome_id", "Genome")
ANNOTATION_FILE_COLUMN = "annotation_file"


def _parse_ko_list_arg(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def _print_summary(result: dict, output_path: str | None) -> None:
    parts = [
        f"unique_ko_count={result.get('unique_ko_count', '?')}",
        f"modules_with_any_hit={result.get('modules_with_any_hit', '?')}",
        f"modules_above_threshold={result.get('modules_above_threshold', len(result.get('modules', [])))}",
        f"completeness_unit={result.get('completeness_unit', 'percent')}",
    ]
    if result.get("cache_hit"):
        parts.append("cache_hit=true")
    if output_path:
        parts.insert(0, f"output={output_path}")
    print(" ".join(parts), file=sys.stderr)


def _run_analyze(args: argparse.Namespace) -> int:
    ko_items = _parse_ko_list_arg(args.ko_list)
    if not args.annotation_file and not ko_items:
        print(
            "error: provide --annotation-file and/or --ko-list with KO content",
            file=sys.stderr,
        )
        return 1

    try:
        kos, ko_meta = resolve_kos(
            ko_list=ko_items or None,
            annotation_file=args.annotation_file or None,
            kegg_column=args.kegg_column,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    analyzer = KeggModuleAnalyzer(KeggAPI())
    ko_sorted = sorted(kos)

    if args.module_ids:
        ids = [m.strip() for m in args.module_ids.split(",") if m.strip()]
        result = analyzer.check_specific_modules(ids, ko_sorted)
    elif args.org:
        result = analyzer.check_organism_modules(
            args.org, ko_sorted, min_completeness=args.min_completeness
        )
    else:
        result = analyzer.check_all_modules(
            ko_sorted, min_completeness=args.min_completeness
        )
    result = {**ko_meta, **result}

    output_path = args.output
    if not output_path:
        default_path = default_tsv_output_path(args.annotation_file or None)
        output_path = str(default_path) if default_path else None

    if output_path:
        path = write_module_completeness_tsv(result["modules"], output_path)
        _print_summary(result, str(path))
    else:
        _print_summary(result, None)
        print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


def _annotation_path_for_genome(
    genome_id: str, downloads_root: str | Path
) -> Path:
    root = Path(downloads_root)
    return root / genome_id / "genome" / f"{genome_id}_gene_annotations.tsv"


def _manifest_data_row_count(manifest_path: Path) -> int:
    with open(manifest_path, "rb") as f:
        return max(0, sum(1 for _ in f) - 1)


def _detect_manifest_columns(fieldnames: list[str] | None) -> tuple[str | None, str | None]:
    if not fieldnames:
        return None, None
    lower_map = {name.lower(): name for name in fieldnames}
    for col in GENOME_ID_COLUMNS:
        if col in fieldnames:
            return col, None
        if col.lower() in lower_map:
            return lower_map[col.lower()], None
    if ANNOTATION_FILE_COLUMN in fieldnames:
        return None, ANNOTATION_FILE_COLUMN
    if "annotation_file" in lower_map:
        return None, lower_map["annotation_file"]
    return fieldnames[0], None


def _load_manifest_rows(
    manifest_path: Path, limit: int | None
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with open(manifest_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        genome_col, annotation_col = _detect_manifest_columns(reader.fieldnames)
        if genome_col is None and annotation_col is None:
            raise ValueError(
                f"Manifest {manifest_path} has no recognizable columns "
                f"(expected genome_id, Genome, annotation_file, or a header row)"
            )
        for row in reader:
            if not any((v or "").strip() for v in row.values()):
                continue
            genome_id = (row.get(genome_col or "", "") or "").strip()
            annotation_file = (row.get(annotation_col or "", "") or "").strip()
            rows.append(
                {
                    "genome_id": genome_id,
                    "annotation_file": annotation_file,
                }
            )
            if limit is not None and len(rows) >= limit:
                break
    return rows


def _resolve_batch_job(
    row: dict[str, str],
    downloads_root: Path,
    output_dir: Path | None,
) -> tuple[str, Path, Path]:
    genome_id = row["genome_id"]
    annotation_file = row.get("annotation_file", "").strip()
    if annotation_file:
        ann_path = Path(annotation_file)
        if not genome_id:
            genome_id = ann_path.name.split("_")[0] if "_" in ann_path.name else ann_path.stem
    else:
        if not genome_id:
            raise ValueError("Manifest row missing genome_id and annotation_file")
        ann_path = _annotation_path_for_genome(genome_id, downloads_root)
    if output_dir is not None:
        out_path = output_dir / f"{genome_id}_module_completeness.tsv"
    else:
        out_path = default_tsv_output_path(ann_path)
        if out_path is None:
            raise ValueError(f"Cannot derive output path for {ann_path}")
    return genome_id, ann_path, out_path


def _batch_job_payload(job: dict) -> dict:
    return {
        "genome_id": job["genome_id"],
        "annotation_file": str(job["annotation_file"]),
        "output_path": str(job["output_path"]),
        "kegg_column": job["kegg_column"],
        "min_completeness": job["min_completeness"],
        "skip_existing": job["skip_existing"],
    }


def _run_batch_worker(payload: dict) -> dict:
    genome_id = payload["genome_id"]
    output_path = Path(payload["output_path"])
    if payload["skip_existing"] and output_path.is_file() and output_path.stat().st_size > 0:
        return {
            "genome_id": genome_id,
            "status": "skip",
            "output_path": str(output_path),
        }
    try:
        kos, ko_meta = resolve_kos(
            annotation_file=payload["annotation_file"],
            kegg_column=payload["kegg_column"],
        )
        analyzer = KeggModuleAnalyzer(KeggAPI())
        result = analyzer.check_all_modules(
            sorted(kos), min_completeness=payload["min_completeness"]
        )
        result = {**ko_meta, **result}
        write_module_completeness_tsv(result["modules"], output_path)
        return {
            "genome_id": genome_id,
            "status": "ok",
            "output_path": str(output_path),
            "unique_ko_count": result.get("unique_ko_count"),
            "modules_with_any_hit": result.get("modules_with_any_hit"),
            "modules_above_threshold": result.get("modules_above_threshold"),
            "cache_hit": result.get("cache_hit", False),
        }
    except Exception as exc:
        return {
            "genome_id": genome_id,
            "status": "fail",
            "error": str(exc),
        }


def _run_batch(args: argparse.Namespace) -> int:
    manifest_path = Path(args.manifest)
    if not manifest_path.is_file():
        print(f"error: manifest not found: {manifest_path}", file=sys.stderr)
        return 1

    data_rows = _manifest_data_row_count(manifest_path)
    if data_rows > MANIFEST_LARGE_ROW_THRESHOLD and args.limit is None:
        print(
            f"error: manifest has {data_rows} data rows (>{MANIFEST_LARGE_ROW_THRESHOLD}). "
            "Filter the manifest or pass --limit N to process a subset.",
            file=sys.stderr,
        )
        return 1

    try:
        rows = _load_manifest_rows(manifest_path, args.limit)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if not rows:
        print("error: manifest has no data rows", file=sys.stderr)
        return 1

    downloads_root = Path(args.downloads_root)
    output_dir = Path(args.output_dir) if args.output_dir else None
    jobs: list[dict] = []
    for row in rows:
        try:
            genome_id, ann_path, out_path = _resolve_batch_job(
                row, downloads_root, output_dir
            )
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        if not ann_path.is_file():
            print(
                f"error: annotation file missing for {genome_id}: {ann_path}",
                file=sys.stderr,
            )
            return 1
        jobs.append(
            {
                "genome_id": genome_id,
                "annotation_file": ann_path,
                "output_path": out_path,
                "kegg_column": args.kegg_column,
                "min_completeness": args.min_completeness,
                "skip_existing": args.skip_existing,
            }
        )

    counts = {"ok": 0, "fail": 0, "skip": 0}
    max_workers = max(1, int(args.jobs))

    def _handle_result(item: dict) -> None:
        status = item.get("status", "fail")
        counts[status] = counts.get(status, 0) + 1
        gid = item.get("genome_id", "?")
        if status == "ok":
            print(
                f"{gid}: ok unique_ko_count={item.get('unique_ko_count')} "
                f"modules_with_any_hit={item.get('modules_with_any_hit')} "
                f"modules_above_threshold={item.get('modules_above_threshold')} "
                f"output={item.get('output_path')}",
                file=sys.stderr,
            )
        elif status == "skip":
            print(f"{gid}: skip existing {item.get('output_path')}", file=sys.stderr)
        else:
            print(f"{gid}: fail {item.get('error', 'unknown')}", file=sys.stderr)

    payloads = [_batch_job_payload(job) for job in jobs]
    if max_workers == 1:
        for payload in payloads:
            _handle_result(_run_batch_worker(payload))
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            futures = [pool.submit(_run_batch_worker, p) for p in payloads]
            for future in as_completed(futures):
                _handle_result(future.result())

    print(
        f"batch done: ok={counts['ok']} fail={counts['fail']} skip={counts['skip']} "
        f"total={len(jobs)}",
        file=sys.stderr,
    )
    return 1 if counts["fail"] else 0


def _run_cog_lookup(args: argparse.Namespace) -> int:
    try:
        catalog = get_cog_catalog(args.cog_csv or None)
        ids = _parse_ko_list_arg(args.cog_ids)
        if not ids:
            print("error: provide --cog-ids", file=sys.stderr)
            return 1
        print(json.dumps(catalog.lookup_many(ids), indent=2, ensure_ascii=False))
        return 0
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _run_cog_search(args: argparse.Namespace) -> int:
    try:
        catalog = get_cog_catalog(args.cog_csv or None)
        hits = catalog.search(args.query, limit=args.limit, cat=args.cat or "")
        print(
            json.dumps(
                {"query": args.query, "count": len(hits), "results": hits},
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="kegg_mcp")
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser(
        "analyze",
        help="Check KEGG module completeness from annotation file and/or KO list",
    )
    analyze.add_argument(
        "--annotation-file",
        default="",
        help="Gene annotation TSV with a KEGG column",
    )
    analyze.add_argument("--kegg-column", default="KEGG", help="KEGG column name")
    analyze.add_argument(
        "--ko-list",
        default="",
        help="Comma-separated KO strings (merged with annotation file)",
    )
    analyze.add_argument(
        "--min-completeness",
        type=float,
        default=0.0,
        help="Minimum completeness threshold (0-100 percent)",
    )
    analyze.add_argument(
        "--module-ids",
        default="",
        help="Comma-separated module IDs (e.g. M00001,M00002)",
    )
    analyze.add_argument(
        "--org",
        default="",
        help="KEGG organism code to restrict module set",
    )
    analyze.add_argument(
        "--output",
        default="",
        help="Output TSV path (default: {genome_id}_module_completeness.tsv beside annotation file)",
    )
    analyze.set_defaults(func=_run_analyze)

    batch = subparsers.add_parser(
        "batch",
        help="Batch module completeness from a manifest TSV",
    )
    batch.add_argument(
        "--manifest",
        required=True,
        help="TSV with genome_id/Genome or annotation_file column",
    )
    batch.add_argument(
        "--downloads-root",
        default="downloads",
        help="Root directory for {genome_id}/genome/*_gene_annotations.tsv",
    )
    batch.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel worker processes (default: 1)",
    )
    batch.add_argument("--kegg-column", default="KEGG", help="KEGG column name")
    batch.add_argument(
        "--min-completeness",
        type=float,
        default=0.0,
        help="Minimum completeness threshold (0-100 percent)",
    )
    batch.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip when output TSV exists and is non-empty",
    )
    batch.add_argument(
        "--output-dir",
        default="",
        help="Write all outputs to this directory (default: beside annotation file)",
    )
    batch.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process at most N manifest rows (required for large manifests)",
    )
    batch.set_defaults(func=_run_batch)

    cog = subparsers.add_parser("cog", help="Query local COG.csv catalog")
    cog_sub = cog.add_subparsers(dest="cog_command", required=True)

    cog_lookup = cog_sub.add_parser("lookup", help="Look up COG IDs")
    cog_lookup.add_argument(
        "--cog-ids",
        required=True,
        help="Comma-separated COG IDs (e.g. COG1387,COG:COG1564)",
    )
    cog_lookup.add_argument(
        "--cog-csv",
        default="",
        help="Path to COG.csv (default: project COG.csv or KEGG_MCP_COG_CSV)",
    )
    cog_lookup.set_defaults(func=_run_cog_lookup)

    cog_search = cog_sub.add_parser("search", help="Search COG annotations")
    cog_search.add_argument("query", help="Keyword in annotation, gene, pathway, or ID")
    cog_search.add_argument("--limit", type=int, default=20)
    cog_search.add_argument("--cat", default="", help="COG category letter filter")
    cog_search.add_argument("--cog-csv", default="")
    cog_search.set_defaults(func=_run_cog_search)

    args = parser.parse_args(argv)
    return args.func(args)
