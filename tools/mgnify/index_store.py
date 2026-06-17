"""Read release_manifest + species_index JSONL; search with match ranking."""

from __future__ import annotations

import csv
import json
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

from .download_species_annotations import normalize_species_name

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = REPO_ROOT / "data" / "mgnify"
HUMAN_GUT_DEFAULT_RELEASE = "v2.0"
_MGYG_RE = re.compile(r"^MGYG[0-9]+$", re.IGNORECASE)

MatchType = Literal["exact", "genus", "fuzzy", "id_resolve"]
_MATCH_RANK = {"exact": 0, "id_resolve": 1, "genus": 2, "fuzzy": 3}


def repo_root() -> Path:
    return REPO_ROOT


def default_data_dir() -> Path:
    env = os.environ.get("MGNIFY_DATA_DIR", "").strip()
    if env:
        return Path(env)
    return DEFAULT_DATA_DIR


def release_sort_key(release: str) -> tuple[int, ...]:
    text = release.lstrip("vV")
    parts: list[int] = []
    for segment in text.split("."):
        try:
            parts.append(int(segment))
        except ValueError:
            parts.append(0)
    return tuple(parts)


def resolve_release(
    manifest: list[dict[str, Any]],
    biome: str,
    release_spec: str = "latest",
) -> str:
    """Resolve release: latest → is_latest, else max version; human-gut fallback v2.0."""
    entries = [e for e in manifest if e.get("biome") == biome]
    if release_spec != "latest":
        return release_spec

    if not entries:
        if biome == "human-gut":
            return HUMAN_GUT_DEFAULT_RELEASE
        raise ValueError(f"release_manifest 中无 biome={biome} 条目")

    latest_flagged = [e for e in entries if e.get("is_latest")]
    pool = latest_flagged if latest_flagged else entries
    best = max(pool, key=lambda e: release_sort_key(str(e.get("release", ""))))
    return str(best["release"])


def species_index_filename(biome: str, release: str) -> str:
    return f"{biome}_{release}_species_index.jsonl"


class IndexStore:
    def __init__(self, data_dir: Path | None = None, repo: Path | None = None) -> None:
        self.repo_root = repo or repo_root()
        self.data_dir = Path(data_dir) if data_dir else default_data_dir()
        self._manifest: list[dict[str, Any]] | None = None

    @property
    def manifest_path(self) -> Path:
        return self.data_dir / "release_manifest.json"

    def load_manifest(self) -> list[dict[str, Any]]:
        if self._manifest is not None:
            return self._manifest
        path = self.manifest_path
        if not path.is_file():
            raise FileNotFoundError(f"release_manifest 不存在: {path}")
        with path.open(encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            self._manifest = [data]
        elif isinstance(data, list):
            self._manifest = data
        else:
            raise ValueError(f"无效的 release_manifest: {path}")
        return self._manifest

    def resolve_release(self, biome: str, release_spec: str = "latest") -> str:
        return resolve_release(self.load_manifest(), biome, release_spec)

    def species_index_path(self, biome: str, release: str) -> Path:
        return self.data_dir / species_index_filename(biome, release)

    def metadata_path_for(self, biome: str, release: str) -> Path | None:
        manifest = self.load_manifest()
        for entry in manifest:
            if entry.get("biome") == biome and entry.get("release") == release:
                raw = entry.get("metadata_local_path")
                if not raw:
                    return None
                path = Path(str(raw))
                if not path.is_absolute():
                    path = self.repo_root / path
                return path if path.is_file() else None
        return None

    def iter_species_entries(self, biome: str, release: str):
        path = self.species_index_path(biome, release)
        if not path.is_file():
            raise FileNotFoundError(f"species_index 不存在: {path}")
        with path.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                entry = json.loads(line)
                if entry.get("biome") == biome and entry.get("release") == release:
                    yield entry

    def get_species_entry(
        self, species_rep: str, biome: str, release: str
    ) -> dict[str, Any] | None:
        for entry in self.iter_species_entries(biome, release):
            if entry.get("species_rep") == species_rep:
                return entry
        return None

    def resolve_species_rep_from_mgyg(
        self, mgyg_id: str, biome: str, release: str
    ) -> str | None:
        mgyg_id = mgyg_id.upper()
        hit = self.get_species_entry(mgyg_id, biome, release)
        if hit:
            return mgyg_id
        meta_path = self.metadata_path_for(biome, release)
        if meta_path is None:
            return None
        with meta_path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("Genome", "") == mgyg_id or row.get("Species_rep", "") == mgyg_id:
                    return row.get("Species_rep", "").strip() or None
        return None

    @staticmethod
    def _classify_match(entry: dict[str, Any], query_norm: str, query_raw: str) -> MatchType | None:
        name_norm = entry.get("species_name_norm", "")
        lineage = entry.get("lineage", "")

        if name_norm == query_norm:
            return "exact"
        if query_norm and (
            name_norm.startswith(query_norm + " ") or name_norm.startswith(query_norm + "_")
        ):
            return "genus"
        if query_raw.casefold() in lineage.casefold():
            return "fuzzy"
        return None

    @staticmethod
    def _sort_key(entry: dict[str, Any], match_type: MatchType) -> tuple:
        genome_type = entry.get("genome_type", "")
        isolate_rank = 0 if genome_type == "Isolate" else 1
        return (
            _MATCH_RANK[match_type],
            -float(entry.get("completeness") or 0),
            float(entry.get("contamination") or 0),
            isolate_rank,
        )

    @staticmethod
    def item_from_entry(entry: dict[str, Any], match_type: MatchType) -> dict[str, Any]:
        return {
            "species_rep": entry["species_rep"],
            "species_name": entry.get("species_name", ""),
            "match_type": match_type,
            "catalogue_url": entry.get("catalogue_url", ""),
            "cache_dir": entry.get("cache_dir", ""),
            "status": entry.get("status", "missing"),
            "completeness": entry.get("completeness"),
            "contamination": entry.get("contamination"),
            "genome_count": entry.get("genome_count"),
        }

    def search_species(
        self,
        biome: str,
        species_query: str,
        release: str = "latest",
        limit: int = 10,
    ) -> dict[str, Any]:
        release_resolved = self.resolve_release(biome, release)
        query_raw = species_query.strip()
        query_norm = normalize_species_name(query_raw)

        if _MGYG_RE.match(query_raw):
            rep = self.resolve_species_rep_from_mgyg(query_raw.upper(), biome, release_resolved)
            if rep is None:
                return {
                    "items": [],
                    "total": 0,
                    "limit": limit,
                    "release": release_resolved,
                }
            entry = self.get_species_entry(rep, biome, release_resolved)
            if entry is None:
                return {
                    "items": [],
                    "total": 0,
                    "limit": limit,
                    "release": release_resolved,
                }
            items = [self.item_from_entry(entry, "id_resolve")]
            return {
                "items": items[:limit],
                "total": len(items),
                "limit": limit,
                "release": release_resolved,
            }

        scored: list[tuple[MatchType, dict[str, Any]]] = []
        for entry in self.iter_species_entries(biome, release_resolved):
            match_type = self._classify_match(entry, query_norm, query_raw)
            if match_type:
                scored.append((match_type, entry))

        scored.sort(key=lambda pair: self._sort_key(pair[1], pair[0]))
        items = [self.item_from_entry(e, mt) for mt, e in scored[:limit]]
        return {
            "items": items,
            "total": len(scored),
            "limit": limit,
            "release": release_resolved,
        }

    def releases_summary(self) -> list[dict[str, Any]]:
        manifest = self.load_manifest()
        return [
            {
                "manifest_id": e.get("manifest_id"),
                "biome": e.get("biome"),
                "release": e.get("release"),
                "is_latest": e.get("is_latest", False),
                "metadata_file_url": e.get("metadata_file_url"),
                "species_index_file": species_index_filename(
                    str(e.get("biome", "")), str(e.get("release", ""))
                ),
            }
            for e in manifest
        ]


@lru_cache(maxsize=1)
def get_default_store() -> IndexStore:
    return IndexStore()
