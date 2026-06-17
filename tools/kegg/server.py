"""FastMCP server for KEGG APIs and module completeness."""

from __future__ import annotations

import json

from mcp.server.fastmcp import FastMCP

from .cog_store import entry_as_json, get_cog_catalog, meta_as_json
from .kegg_api import KeggAPI
from .ko_input import resolve_kos
from .module_completeness import (
    KeggModuleAnalyzer,
    default_tsv_output_path,
    write_module_completeness_tsv,
)
from .allowlist_sources import (
    BLACKLIST_VERSION,
    blacklist_summary,
    all_blacklist_entries,
    load_allowlist,
    load_blacklist,
    write_blacklist_files,
)

mcp = FastMCP("kegg_mcp", stateless_http=True)
kegg_api = KeggAPI()
module_analyzer = KeggModuleAnalyzer(kegg_api)


def _parse_csv_arg(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


@mcp.tool()
async def kegg_info(database: str):
    """Get KEGG database statistics via info operation."""
    try:
        text = kegg_api.kegg_info(database)
        return [{"result": text}]
    except Exception as e:
        return [{"error": f"An error occurred while querying KEGG info: {str(e)}"}]


@mcp.tool()
async def kegg_find(database: str, query, option: str | None = None):
    """Find entries by keyword/formula/mass in KEGG."""
    try:
        result = kegg_api.kegg_find(database, query, option)
    except Exception as e:
        return {"error": f"An error occurred while querying KEGG find: {str(e)}"}
    return {"result": result}


@mcp.tool()
async def kegg_list(database, org: str | None = None):
    """List KEGG entries in a database."""
    try:
        result = kegg_api.kegg_list(database, org)
    except Exception as e:
        return {"error": f"An error occurred while querying KEGG list: {str(e)}"}
    return {"result": result}


@mcp.tool()
async def kegg_get(dbentries, option: str | None = None):
    """Get KEGG flat-file/raw content by entry ID(s)."""
    try:
        result = kegg_api.kegg_get(dbentries, option)
    except Exception as e:
        return {"error": f"An error occurred while querying KEGG get: {str(e)}"}
    return {"result": result}


@mcp.tool()
async def kegg_conv(target_db: str, source_db, option: str | None = None):
    """Convert identifiers between KEGG and external DBs."""
    try:
        result = kegg_api.kegg_conv(target_db, source_db, option)
    except Exception as e:
        return {"error": f"An error occurred while querying KEGG conv: {str(e)}"}
    return {"result": result}


@mcp.tool()
async def kegg_link(target_db: str, source_db, option: str | None = None):
    """Link entries across KEGG databases."""
    try:
        result = kegg_api.kegg_link(target_db, source_db, option)
    except Exception as e:
        return {"error": f"An error occurred while querying KEGG link: {str(e)}"}
    return {"result": result}


@mcp.tool()
async def kegg_module_completeness(
    ko_list: str = "",
    annotation_file: str = "",
    kegg_column: str = "KEGG",
    module_ids: str = "",
    org: str = "",
    min_completeness: float = 0.0,
    output_format: str = "json",
    output_path: str = "",
):
    """
    Check metabolic module completeness against KO identifiers from an annotation
    file and/or a comma-separated ko_list.

    Provide annotation_file and/or ko_list (at least one with valid KOs).
    When module_ids and org are both empty, scans all modules with completeness > 0
    (mode=all). completeness values are percentages (0-100); see completeness_unit.
    Response metadata includes unique_ko_count, modules_with_any_hit,
    modules_above_threshold, tool, and tool_version.

    Modes:
    1. All modules (default): module_ids and org empty -> check_all_modules
    2. Specific modules: module_ids comma-separated (e.g. "M00001,M00002")
    3. Organism subset: org KEGG code (e.g. "hsa") when module_ids empty
    """
    try:
        ko_items = _parse_csv_arg(ko_list)
        if not (annotation_file.strip() or ko_items):
            return {
                "error": (
                    "At least one of annotation_file or ko_list must be provided"
                )
            }
        kos, ko_meta = resolve_kos(
            ko_list=ko_items or None,
            annotation_file=annotation_file.strip() or None,
            kegg_column=kegg_column,
        )
        ko_sorted = sorted(kos)

        ids = _parse_csv_arg(module_ids)
        if ids:
            result = module_analyzer.check_specific_modules(ids, ko_sorted)
        elif org.strip():
            result = module_analyzer.check_organism_modules(
                org.strip(), ko_sorted, min_completeness
            )
        else:
            result = module_analyzer.check_all_modules(ko_sorted, min_completeness)
        result = {**ko_meta, **result}

        if output_format.lower() == "tsv":
            tsv_path = output_path.strip()
            if not tsv_path:
                default_path = default_tsv_output_path(annotation_file or None)
                if default_path is None:
                    return {
                        "error": (
                            "output_path required for TSV when annotation_file "
                            "is not provided"
                        )
                    }
                tsv_path = str(default_path)
            written = write_module_completeness_tsv(result["modules"], tsv_path)
            result = {**result, "output_tsv": str(written)}
    except Exception as e:
        return {"error": f"An error occurred while checking module completeness: {str(e)}"}
    return result


@mcp.resource(
    "resource://cog/catalog/meta",
    name="cog_catalog_meta",
    description="COG catalog metadata (KEGG_MCP_COG_CSV, kegg_mcp/data/COG.csv, or project COG.csv).",
    mime_type="application/json",
)
def cog_catalog_meta() -> str:
    """COG catalog index metadata (resource://cog/catalog/meta)."""
    try:
        return meta_as_json()
    except Exception as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)


@mcp.resource(
    "resource://cog/entry/{cog_id}",
    name="cog_entry",
    description="Single COG term by ID (requires COG catalog file; see cog/catalog/meta).",
    mime_type="application/json",
)
def cog_entry_resource(cog_id: str) -> str:
    """One COG record by ID, e.g. resource://cog/entry/COG1387."""
    try:
        catalog = get_cog_catalog()
        return entry_as_json(catalog.lookup(cog_id), cog_id)
    except Exception as e:
        return json.dumps({"error": str(e), "cog_id": cog_id}, ensure_ascii=False)


@mcp.tool()
async def cog_lookup(cog_ids: str):
    """
    Look up COG term descriptions from the configured COG catalog file.

    Accepts comma-separated IDs in forms like COG1387, cog:COG1387, or COG:COG1387.
    Returns a map of normalized COG ID -> entry dict (or null if not found).
    """
    try:
        ids = _parse_csv_arg(cog_ids)
        if not ids:
            return {"error": "cog_ids is required (comma-separated COG identifiers)"}
        catalog = get_cog_catalog()
        return {
            "source_file": str(catalog.csv_path),
            "results": catalog.lookup_many(ids),
        }
    except Exception as e:
        return {"error": f"COG lookup failed: {e}"}


@mcp.tool()
async def cog_search(query: str, limit: int = 20, cat: str = ""):
    """
    Keyword search in the COG catalog (annotation, gene symbol, pathway, COG ID).

    Optional cat filters by COG category letter(s), e.g. E or GH.
    """
    try:
        catalog = get_cog_catalog()
        hits = catalog.search(query, limit=limit, cat=cat)
        return {
            "source_file": str(catalog.csv_path),
            "query": query,
            "count": len(hits),
            "results": hits,
        }
    except Exception as e:
        return {"error": f"COG search failed: {e}"}


@mcp.tool()
async def kegg_module_detail(module_id: str):
    """
    Retrieve detailed metadata and definition for a KEGG metabolic module.

    Query example: {"module_id": "M00001"}
    """
    try:
        result = module_analyzer.fetch_module_definition(module_id)
    except Exception as e:
        return {"error": f"An error occurred while fetching module detail: {str(e)}"}
    return result


# ---------------------------------------------------------------------------
# Allowlist / Blacklist Resources
# ---------------------------------------------------------------------------

@mcp.resource(
    "resource://allowlist/meta",
    name="allowlist_meta",
    description="Prokaryote allowlist metadata (versions, domains, statistics).",
    mime_type="application/json",
)
def allowlist_meta() -> str:
    """Allowlist metadata index (resource://allowlist/meta)."""
    try:
        from pathlib import Path
        data_dir = Path(__file__).resolve().parent / "data"

        versions = ["2.0", "0.1.2", "0.1.1", "0.1"]
        meta_list = []

        for ver in versions:
            split_meta_path = data_dir / f"prokaryote_allowlist_split_meta_v{ver}.json"
            legacy_meta_path = data_dir / f"prokaryote_allowlist_meta_v{ver}.json"

            if split_meta_path.exists():
                content = json.loads(split_meta_path.read_text(encoding="utf-8"))
                content["type"] = "split"
                meta_list.append(content)
            elif legacy_meta_path.exists():
                content = json.loads(legacy_meta_path.read_text(encoding="utf-8"))
                content["type"] = "legacy"
                meta_list.append(content)

        return json.dumps({
            "versions": [m["version"] for m in meta_list],
            "allowlists": meta_list,
        }, indent=2, ensure_ascii=False)
    except Exception as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)


@mcp.resource(
    "resource://allowlist/{domain}/{kind}/{version}",
    name="allowlist_entries",
    description="Prokaryote allowlist entries (domain: bacteria/archaea/prok, kind: module/pathway, version: 2.0/0.1.2/0.1.1/0.1).",
    mime_type="application/json",
)
def allowlist_entries(domain: str, kind: str, version: str) -> str:
    """Allowlist entries by domain, kind, and version (resource://allowlist/{domain}/{kind}/{version})."""
    try:
        from pathlib import Path
        data_dir = Path(__file__).resolve().parent / "data"

        allowed_domains = {"bacteria", "archaea", "prok"}
        allowed_kinds = {"module", "pathway"}
        allowed_versions = {"2.0", "0.1.2", "0.1.1", "0.1"}

        if domain not in allowed_domains:
            return json.dumps({"error": f"domain must be one of {allowed_domains}"}, ensure_ascii=False)
        if kind not in allowed_kinds:
            return json.dumps({"error": f"kind must be one of {allowed_kinds}"}, ensure_ascii=False)
        if version not in allowed_versions:
            return json.dumps({"error": f"version must be one of {allowed_versions}"}, ensure_ascii=False)

        ids = load_allowlist(data_dir, domain, kind, version)

        # Load full metadata from TSV
        import csv
        metadata: dict[str, dict] = {}

        if kind == "module":
            tsv_path = data_dir / f"kegg_module_allowlist_{domain}_v{version}.tsv"
            id_col = "module_id"
            name_col = "module_name"
        else:
            tsv_path = data_dir / f"prokaryote_ko_pathway_allowlist_{domain}_v{version}.tsv"
            id_col = "pathway_id"
            name_col = "pathway_name"

        if tsv_path.exists():
            with tsv_path.open(encoding="utf-8") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    metadata[row[id_col]] = row

        entries = []
        for entry_id in sorted(ids):
            meta = metadata.get(entry_id, {})
            entries.append({
                "id": entry_id,
                "name": meta.get(name_col, ""),
                **{k: v for k, v in meta.items() if k not in (id_col, name_col)},
            })

        return json.dumps({
            "domain": domain,
            "kind": kind,
            "version": version,
            "count": len(entries),
            "entries": entries,
        }, indent=2, ensure_ascii=False)
    except Exception as e:
        return json.dumps({"error": str(e), "domain": domain, "kind": kind, "version": version}, ensure_ascii=False)


@mcp.resource(
    "resource://blacklist/meta",
    name="blacklist_meta",
    description="Prokaryote blacklist metadata (statistics, categories, rules).",
    mime_type="application/json",
)
def blacklist_meta() -> str:
    """Blacklist metadata and summary (resource://blacklist/meta)."""
    try:
        summary = blacklist_summary()
        return json.dumps(summary, indent=2, ensure_ascii=False)
    except Exception as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)


@mcp.resource(
    "resource://blacklist/{target}/{version}",
    name="blacklist_entries",
    description="Blacklist entries for module or pathway (target: module/pathway, version: 1.0).",
    mime_type="application/json",
)
def blacklist_entries(target: str, version: str = BLACKLIST_VERSION) -> str:
    """Blacklist entries by target and version (resource://blacklist/{target}/{version})."""
    try:
        from dataclasses import asdict

        allowed_targets = {"module", "pathway"}
        if target not in allowed_targets:
            return json.dumps({"error": f"target must be one of {allowed_targets}"}, ensure_ascii=False)

        entries = all_blacklist_entries()
        filtered = [e for e in entries if e.applies_to == "both" or e.applies_to == target]

        # Add ID/prefix rules for pathway
        if target == "pathway":
            from .allowlist_sources import PATHWAY_EXCLUDE_IDS, PATHWAY_EXCLUDE_PREFIXES
            for ko_id in sorted(PATHWAY_EXCLUDE_IDS):
                filtered.append({"keyword": ko_id, "category": "global_map", "reason": "全局代谢图，非独立通路", "rule_type": "id"})
            for prefix in PATHWAY_EXCLUDE_PREFIXES:
                filtered.append({"keyword": prefix, "category": "pathway_block", "reason": "KEGG通路分类块（细胞过程/生物体系统/人类疾病）真核为主", "rule_type": "prefix"})

        return json.dumps({
            "target": target,
            "version": version,
            "count": len(filtered),
            "entries": [asdict(e) if hasattr(e, "__dataclass_fields__") else e for e in filtered],
        }, indent=2, ensure_ascii=False)
    except Exception as e:
        return json.dumps({"error": str(e), "target": target, "version": version}, ensure_ascii=False)


if __name__ == "__main__":
    mcp.run(transport="stdio")
