"""KEGG module completeness analyzer backed by kegg-pathways-completeness CLI."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from collections.abc import Callable
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .kegg_api import KeggAPI
from .ko_input import normalize_kos

TOOL_NAME = "kegg-pathways-completeness"


def _tool_version() -> str:
    try:
        return version(TOOL_NAME)
    except PackageNotFoundError:
        return "unknown"


def _enrich_result_metadata(
    result: dict,
    kos: set[str],
    min_completeness: float = 0.0,
    *,
    cache_hit: bool = False,
) -> dict:
    modules = result.get("modules", [])
    if "modules_with_any_hit" not in result:
        result["modules_with_any_hit"] = sum(
            1 for item in modules if float(item.get("completeness", 0.0)) > 0
        )
    result["unique_ko_count"] = len(kos)
    result["ko_input_count"] = len(kos)
    result["modules_above_threshold"] = len(modules)
    result["completeness_unit"] = "percent"
    result["tool"] = TOOL_NAME
    result["tool_version"] = _tool_version()
    if cache_hit:
        result["cache_hit"] = True
    return result

TSV_COLUMNS = (
    "module_id",
    "module_name",
    "module_class",
    "completeness",
    "matching_ko",
    "missing_ko",
)


class KeggModuleAnalyzer:
    """Analyze KEGG module completeness from KO list."""

    _completeness_cache: dict[frozenset[str], dict[str, dict[str, str]]] = {}

    def __init__(self, kegg_api: KeggAPI):
        self.kegg_api = kegg_api
        self._last_cache_hit = False

    @classmethod
    def clear_cache(cls) -> None:
        cls._completeness_cache.clear()

    def fetch_module_definition(self, module_id: str) -> dict:
        flat_text = self.kegg_api.kegg_get(module_id)
        if not flat_text.strip():
            raise ValueError(f"Failed to fetch module: {module_id}")

        fields: dict[str, str] = {}
        current_key: str | None = None
        for raw_line in flat_text.splitlines():
            if not raw_line.strip():
                continue
            key = raw_line[:12].strip()
            value = raw_line[12:].strip()
            if key:
                current_key = key
                if key == "DEFINITION":
                    fields[key] = (fields.get(key, "") + "\n" + value).strip()
                else:
                    fields[key] = value
            elif current_key:
                joiner = "\n" if current_key == "DEFINITION" else " "
                fields[current_key] = f"{fields.get(current_key, '')}{joiner}{value}".strip()

        definition = fields.get("DEFINITION", "").strip()
        if not definition:
            raise ValueError(f"Missing DEFINITION for module: {module_id}")

        entry = fields.get("ENTRY", module_id).split()[0]
        return {
            "module_id": entry,
            "name": fields.get("NAME", "").strip(),
            "definition": definition,
            "class": fields.get("CLASS", "").strip(),
        }

    def fetch_organism_modules(self, org: str) -> list[dict]:
        listing = self.kegg_api.kegg_list("module", org)
        modules: list[dict] = []
        for line in listing.splitlines():
            first_col = line.split("\t", 1)[0].strip()
            if first_col.startswith("md:"):
                mid = first_col.split(":", 1)[1].strip()
                modules.append({"module_id": mid})
        return modules

    def parse_definition_to_graph(self, definition: str):
        raise NotImplementedError(
            "Graph parsing is delegated to kegg-pathways-completeness."
        )

    def calculate_completeness(self, graph, ko_set: set[str]) -> dict:
        raise NotImplementedError(
            "Completeness is delegated to kegg-pathways-completeness."
        )

    def check_module(self, module_id: str, ko_list: list[str]) -> dict:
        clean_kos = self._normalize_kos(ko_list)
        result_map = self._run_completeness_tool(clean_kos)
        module = self.fetch_module_definition(module_id)
        row = result_map.get(module["module_id"], {})
        completeness = float(row.get("completeness", 0.0))
        result = {
            "module_id": module["module_id"],
            "module_name": module["name"],
            "module_class": module["class"],
            "definition": module["definition"],
            "completeness": completeness,
            "matching_kos": self._split_csv_field(row.get("matching_ko", "")),
            "missing_kos": self._split_csv_field(row.get("missing_ko", "")),
            "optional_missing_kos": [],
            "modules": [
                {
                    "module_id": module["module_id"],
                    "module_name": module["name"],
                    "module_class": module["class"],
                    "completeness": completeness,
                    "matching_kos": self._split_csv_field(row.get("matching_ko", "")),
                    "missing_kos": self._split_csv_field(row.get("missing_ko", "")),
                    "optional_missing_kos": [],
                }
            ],
            "modules_with_any_hit": 1 if completeness > 0 else 0,
        }
        return _enrich_result_metadata(
            result, clean_kos, cache_hit=self._last_cache_hit
        )

    def check_organism_modules(
        self, org: str, ko_list: list[str], min_completeness: float = 0.0
    ) -> dict:
        modules = self.fetch_organism_modules(org)
        module_ids = {m["module_id"] for m in modules}
        clean_kos = self._normalize_kos(ko_list)
        result_map = self._run_completeness_tool(clean_kos)
        results = self._rows_from_result_map(
            result_map,
            lambda module_id, _row: module_id in module_ids,
        )
        results.sort(key=lambda x: x["completeness"], reverse=True)
        modules_with_any_hit = len(results)
        filtered = [item for item in results if item["completeness"] >= min_completeness]
        return _enrich_result_metadata(
            {
                "organism": org,
                "total_modules_checked": len(results),
                "modules_with_any_hit": modules_with_any_hit,
                "modules": filtered,
            },
            clean_kos,
            min_completeness,
            cache_hit=self._last_cache_hit,
        )

    def check_specific_modules(self, module_ids: list[str], ko_list: list[str]) -> dict:
        clean_kos = self._normalize_kos(ko_list)
        wanted = {m.strip() for m in module_ids if m.strip()}
        result_map = self._run_completeness_tool(clean_kos)
        results = self._rows_from_result_map(
            result_map,
            lambda module_id, _row: module_id in wanted,
        )
        found_ids = {item["module_id"] for item in results}
        for module_id in sorted(wanted - found_ids):
            results.append(self._module_entry(module_id, None))
        results.sort(key=lambda x: x["completeness"], reverse=True)
        modules_with_any_hit = sum(
            1 for item in results if float(item.get("completeness", 0.0)) > 0
        )
        return _enrich_result_metadata(
            {
                "total_modules_checked": len(results),
                "modules_with_any_hit": modules_with_any_hit,
                "modules": results,
            },
            clean_kos,
            cache_hit=self._last_cache_hit,
        )

    def check_all_modules(
        self, ko_list: list[str], min_completeness: float = 0.0
    ) -> dict:
        clean_kos = self._normalize_kos(ko_list)
        result_map = self._run_completeness_tool(clean_kos)
        results = self._rows_from_result_map(
            result_map,
            lambda _module_id, row: float(row.get("completeness", 0.0)) > 0,
        )
        modules_with_any_hit = len(results)
        filtered = [item for item in results if item["completeness"] >= min_completeness]
        filtered.sort(key=lambda x: x["completeness"], reverse=True)
        return _enrich_result_metadata(
            {
                "modules_with_any_hit": modules_with_any_hit,
                "modules": filtered,
            },
            clean_kos,
            min_completeness,
            cache_hit=self._last_cache_hit,
        )

    def _rows_from_result_map(
        self,
        result_map: dict[str, dict[str, str]],
        filter_fn: Callable[[str, dict[str, str]], bool],
    ) -> list[dict]:
        results: list[dict] = []
        for module_id, row in result_map.items():
            if not filter_fn(module_id, row):
                continue
            results.append(self._module_entry(module_id, row))
        return results

    def _module_entry(
        self, module_id: str, row: dict[str, str] | None
    ) -> dict:
        if row is None:
            return {
                "module_id": module_id,
                "module_name": "",
                "module_class": "",
                "completeness": 0.0,
                "matching_kos": [],
                "missing_kos": [],
                "optional_missing_kos": [],
            }
        return {
            "module_id": module_id,
            "module_name": row.get("pathway_name", ""),
            "module_class": row.get("pathway_class", ""),
            "completeness": float(row.get("completeness", 0.0)),
            "matching_kos": self._split_csv_field(row.get("matching_ko", "")),
            "missing_kos": self._split_csv_field(row.get("missing_ko", "")),
            "optional_missing_kos": [],
        }

    def _normalize_kos(self, ko_list: list[str]) -> set[str]:
        return normalize_kos(ko_list)

    def _split_csv_field(self, value: str) -> list[str]:
        return sorted({item.strip() for item in value.split(",") if item.strip()})

    def _run_completeness_tool(self, ko_set: set[str]) -> dict[str, dict[str, str]]:
        key = frozenset(ko_set)
        cached = self._completeness_cache.get(key)
        if cached is not None:
            self._last_cache_hit = True
            return cached
        self._last_cache_hit = False
        result_map = self._invoke_completeness_tool(ko_set)
        self._completeness_cache[key] = result_map
        return result_map

    def _invoke_completeness_tool(self, ko_set: set[str]) -> dict[str, dict[str, str]]:
        with tempfile.TemporaryDirectory(prefix="kegg_mcp_") as tmpdir:
            ko_file = f"{tmpdir}/kos.txt"
            outprefix = "kegg_mcp"
            with open(ko_file, "w", encoding="utf-8") as f:
                f.write(",".join(sorted(ko_set)))

            cmd = [
                sys.executable,
                "-m",
                "kegg_pathways_completeness.bin.give_completeness",
                "--input-list",
                ko_file,
                "--list-separator",
                ",",
                "--outdir",
                tmpdir,
                "--outprefix",
                outprefix,
            ]
            run = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if run.returncode != 0:
                stderr = (run.stderr or "").strip()
                stdout = (run.stdout or "").strip()
                raise RuntimeError(
                    "kegg-pathways-completeness failed: "
                    f"{stderr or stdout or 'unknown error'}"
                )

            pathways_file = f"{tmpdir}/{outprefix}_pathways.tsv"
            result_map: dict[str, dict[str, str]] = {}
            with open(pathways_file, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    module_id = (row.get("module_accession") or "").strip()
                    if module_id:
                        result_map[module_id] = row
            return result_map


def write_module_completeness_tsv(modules: list[dict], output_path: str | Path) -> Path:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for item in modules:
            writer.writerow(
                {
                    "module_id": item.get("module_id", ""),
                    "module_name": item.get("module_name", ""),
                    "module_class": item.get("module_class", ""),
                    "completeness": item.get("completeness", 0.0),
                    "matching_ko": ",".join(item.get("matching_kos") or []),
                    "missing_ko": ",".join(item.get("missing_kos") or []),
                }
            )
    return out


def default_tsv_output_path(annotation_file: str | Path | None) -> Path | None:
    if not annotation_file:
        return None
    path = Path(annotation_file)
    genome_id = path.name.split("_")[0] if "_" in path.name else path.stem
    return path.parent / f"{genome_id}_module_completeness.tsv"
