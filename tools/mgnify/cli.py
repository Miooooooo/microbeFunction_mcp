"""CLI: python -m mgnify_mcp search|fetch ..."""

from __future__ import annotations

import argparse
import json
import sys

from .fetch import fetch_annotations
from .index_store import IndexStore


def _run_search(args: argparse.Namespace) -> int:
    store = IndexStore()
    result = store.search_species(
        biome=args.biome,
        species_query=args.query,
        release=args.release,
        limit=args.limit,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


def _run_fetch(args: argparse.Namespace) -> int:
    store = IndexStore()
    roles = [r.strip() for r in args.roles.split(",") if r.strip()] if args.roles else None
    try:
        result = fetch_annotations(
            store,
            species_rep=args.species_rep,
            biome=args.biome,
            release=args.release,
            roles=roles,
            convert_gtf=not args.no_convert_gtf,
            preview_rows=args.preview_rows,
        )
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="mgnify_mcp")
    sub = parser.add_subparsers(dest="command", required=True)

    search = sub.add_parser("search", help="Search species_index")
    search.add_argument("--biome", required=True)
    search.add_argument("--query", required=True, help="Species name or MGYG ID")
    search.add_argument("--release", default="latest")
    search.add_argument("--limit", type=int, default=10)
    search.set_defaults(func=_run_search)

    fetch = sub.add_parser("fetch", help="Fetch annotation files")
    fetch.add_argument("--species-rep", required=True)
    fetch.add_argument("--biome", required=True)
    fetch.add_argument("--release", default="latest")
    fetch.add_argument(
        "--roles",
        default="",
        help="Comma-separated roles (default: gff,eggnog_tsv)",
    )
    fetch.add_argument("--preview-rows", type=int, default=5)
    fetch.add_argument(
        "--no-convert-gtf",
        action="store_true",
        help="Skip GFF→GTF conversion",
    )
    fetch.set_defaults(func=_run_fetch)

    args = parser.parse_args(argv)
    return args.func(args)
