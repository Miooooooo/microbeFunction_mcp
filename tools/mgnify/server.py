"""FastMCP server for MGnify species index search and annotation fetch."""

from __future__ import annotations

import json
from typing import Any

from mcp.server.fastmcp import FastMCP

from . import fetch as fetch_mod
from .index_store import IndexStore, get_default_store

mcp = FastMCP("mgnify_mcp", stateless_http=True)


def _store() -> IndexStore:
    return get_default_store()


@mcp.tool()
async def search_species(
    biome: str,
    species_query: str,
    release: str = "latest",
    limit: int = 10,
) -> dict[str, Any]:
    """
    Search MGnify species_index by species name or MGYG ID.

    release=latest resolves via release_manifest (is_latest, else max version).
    """
    try:
        return _store().search_species(
            biome=biome,
            species_query=species_query,
            release=release,
            limit=limit,
        )
    except Exception as e:
        return {"error": str(e)}


@mcp.tool()
async def fetch_annotations(
    species_rep: str,
    biome: str,
    release: str = "latest",
    roles: list[str] | None = None,
    convert_gtf: bool = True,
    preview_rows: int = 5,
) -> dict[str, Any]:
    """
    Fetch GFF/eggNOG (and optional roles) for a species_rep.

    Uses local cache under downloads/{species_rep}/genome/ when present;
    otherwise downloads from catalogue_url. Returns paths and TSV preview only
    (never full annotation tables).
    """
    try:
        return fetch_mod.fetch_annotations(
            _store(),
            species_rep=species_rep,
            biome=biome,
            release=release,
            roles=roles or ["gff", "eggnog_tsv"],
            convert_gtf=convert_gtf,
            preview_rows=preview_rows,
        )
    except Exception as e:
        return {"error": str(e)}


@mcp.resource(
    "resource://mgnify/releases",
    name="mgnify_releases",
    description="MGnify release_manifest summary (biome, release, is_latest).",
    mime_type="application/json",
)
def mgnify_releases_resource() -> str:
    """Release manifest summary; full species_index lives in local JSONL per biome/release."""
    try:
        store = _store()
        payload = {
            "data_dir": str(store.data_dir),
            "manifest_path": str(store.manifest_path),
            "releases": store.releases_summary(),
            "note": "物种导航见 data/mgnify/{biome}_{release}_species_index.jsonl",
        }
        return json.dumps(payload, ensure_ascii=False, indent=2)
    except Exception as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)


@mcp.prompt()
def system_prompt():
    return """You can search MGnify human-gut (and other biome) species catalogues via
search_species, then fetch annotation file paths and small TSV previews with fetch_annotations.

Do not expect full GFF or eggNOG tables in tool responses—use cache_path for downstream analysis.
Species representative ID is species_rep (not member Genome ID); resolve member IDs via search."""


if __name__ == "__main__":
    mcp.run(transport="stdio")
