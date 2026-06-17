"""Fetch species annotation files from cache or remote FTP."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Literal

from .download_species_annotations import gff_to_gtf
from .list_mgnify_folders import list_files

try:
    from .download_species_annotations import _download_file
except ImportError:
    _download_file = None  # type: ignore

from .index_store import IndexStore

SourceKind = Literal["cache", "remote", "mixed"]
FileStatus = Literal["ready", "missing", "error"]


def _cache_dir_path(repo_root: Path, cache_dir: str) -> Path:
    path = Path(cache_dir)
    if not path.is_absolute():
        path = repo_root / path
    return path


def _basename_candidates(species_rep: str, role: str) -> list[str]:
    mapping: dict[str, list[str]] = {
        "gff": [f"{species_rep}.gff", f"{species_rep}.gff.gz"],
        "gtf": [f"{species_rep}.gtf"],
        "eggnog_tsv": [f"{species_rep}_eggNOG.tsv"],
        "kegg_tsv": [],
        "interpro_tsv": [],
        "gene_annotations_tsv": [f"{species_rep}_gene_annotations.tsv"],
    }
    if role == "kegg_tsv":
        return [f"{species_rep}_KEGG.tsv", f"{species_rep}_kegg.tsv"]
    if role == "interpro_tsv":
        return [
            f"{species_rep}_InterProScan.tsv",
            f"{species_rep}_interpro.tsv",
        ]
    return mapping.get(role, [])


def _find_local_file(cache_path: Path, candidates: list[str]) -> Path | None:
    if not cache_path.is_dir():
        return None
    names = {p.name for p in cache_path.iterdir() if p.is_file()}
    for name in candidates:
        if name in names:
            return cache_path / name
    if candidates and candidates[0].endswith(".gff"):
        gtf = cache_path / candidates[0].replace(".gff", ".gtf")
        if gtf.is_file():
            return None
    return None


def _match_remote_name(remote_files: list[str], candidates: list[str], role: str) -> str | None:
    lower_map = {n.lower(): n for n in remote_files}
    for cand in candidates:
        if cand in remote_files:
            return cand
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    if role == "kegg_tsv":
        for name in remote_files:
            if "kegg" in name.lower() and name.lower().endswith(".tsv"):
                return name
    if role == "interpro_tsv":
        for name in remote_files:
            if "interpro" in name.lower() and name.lower().endswith(".tsv"):
                return name
    return None


def tsv_preview(path: Path, preview_rows: int) -> dict[str, Any]:
    with path.open(encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            columns = next(reader)
        except StopIteration:
            return {"columns": [], "rows": []}
        rows: list[list[str]] = []
        for i, row in enumerate(reader):
            if i >= preview_rows:
                break
            rows.append(row)
    return {"columns": columns, "rows": rows}


def fetch_annotations(
    store: IndexStore,
    *,
    species_rep: str,
    biome: str,
    release: str = "latest",
    roles: list[str] | None = None,
    convert_gtf: bool = True,
    preview_rows: int = 5,
) -> dict[str, Any]:
    roles = roles or ["gff", "eggnog_tsv"]
    release_resolved = store.resolve_release(biome, release)
    entry = store.get_species_entry(species_rep, biome, release_resolved)
    if entry is None:
        raise ValueError(
            f"species_index 中未找到 {species_rep} (biome={biome}, release={release_resolved})"
        )

    catalogue_url = str(entry.get("catalogue_url", "")).rstrip("/") + "/"
    cache_path = _cache_dir_path(store.repo_root, str(entry.get("cache_dir", "")))
    cache_path.mkdir(parents=True, exist_ok=True)

    files_out: list[dict[str, Any]] = []
    preview: dict[str, Any] = {}
    sources: set[str] = set()
    remote_files: list[str] | None = None

    def ensure_remote_list() -> list[str]:
        nonlocal remote_files
        if remote_files is None:
            remote_files = list_files(catalogue_url)
        return remote_files

    for role in roles:
        if role == "gtf":
            continue
        candidates = _basename_candidates(species_rep, role)
        local = _find_local_file(cache_path, candidates)
        status: FileStatus = "missing"

        if local is not None and local.is_file():
            status = "ready"
            sources.add("cache")
        else:
            remote = ensure_remote_list()
            remote_name = _match_remote_name(remote, candidates, role)
            if remote_name is None:
                files_out.append(
                    {
                        "role": role,
                        "cache_path": "",
                        "size_bytes": 0,
                        "status": "missing",
                    }
                )
                continue
            dest = cache_path / remote_name
            if _download_file is None:
                raise RuntimeError("download_species_annotations._download_file 不可用")
            _download_file(catalogue_url + remote_name, dest)
            local = dest
            status = "ready"
            sources.add("remote")

        size_bytes = local.stat().st_size if local and local.is_file() else 0
        rel_path = str(local.relative_to(store.repo_root)).replace("\\", "/")
        files_out.append(
            {
                "role": role,
                "cache_path": rel_path,
                "size_bytes": size_bytes,
                "status": status,
            }
        )

        if preview_rows > 0 and role.endswith("_tsv") and local and local.is_file():
            preview[role] = tsv_preview(local, preview_rows)

        if role == "gff" and convert_gtf and local and local.is_file():
            gtf_path = local.with_suffix(".gtf")
            if not gtf_path.is_file():
                gff_to_gtf(local, gtf_path)
                sources.add("cache")
            files_out.append(
                {
                    "role": "gtf",
                    "cache_path": str(gtf_path.relative_to(store.repo_root)).replace(
                        "\\", "/"
                    ),
                    "size_bytes": gtf_path.stat().st_size,
                    "status": "ready",
                }
            )

    if len(sources) == 1:
        source: SourceKind = next(iter(sources))  # type: ignore[assignment]
    elif len(sources) > 1:
        source = "mixed"
    elif any(f["status"] == "ready" for f in files_out):
        source = "cache"
    else:
        source = "remote"

    result: dict[str, Any] = {
        "species_rep": species_rep,
        "source": source,
        "files": files_out,
        "release": release_resolved,
    }
    if preview:
        result["preview"] = preview
    return result
