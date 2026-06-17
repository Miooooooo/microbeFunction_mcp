"""Tests for mgnify_mcp search and fetch (cache when available)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from tools.mgnify.build_mgnify_index import REFERENCE_SPECIES
from tools.mgnify.fetch import fetch_annotations
from tools.mgnify.index_store import IndexStore, resolve_release

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data" / "mgnify"
BIOME = "human-gut"
RELEASE = "v2.0"
CACHED_SPECIES = "MGYG000000238"
CACHE_DIR = REPO_ROOT / "downloads" / CACHED_SPECIES / "genome"


@pytest.fixture(scope="module")
def store() -> IndexStore:
    if not (DATA_DIR / "release_manifest.json").is_file():
        pytest.skip("release_manifest 未构建")
    return IndexStore(data_dir=DATA_DIR, repo=REPO_ROOT)


def test_resolve_latest_when_is_latest_false(store: IndexStore) -> None:
    manifest = store.load_manifest()
    resolved = resolve_release(manifest, BIOME, "latest")
    assert resolved == RELEASE


@pytest.mark.parametrize("species_name,expected_rep", REFERENCE_SPECIES)
def test_search_reference_species(
    store: IndexStore, species_name: str, expected_rep: str
) -> None:
    result = store.search_species(
        biome=BIOME, species_query=species_name, release=RELEASE, limit=5
    )
    assert result["release"] == RELEASE
    assert result["total"] >= 1
    assert result["items"][0]["species_rep"] == expected_rep
    assert result["items"][0]["match_type"] == "exact"


def test_search_by_species_rep_id(store: IndexStore) -> None:
    result = store.search_species(
        biome=BIOME, species_query=CACHED_SPECIES, release=RELEASE, limit=1
    )
    assert result["items"][0]["species_rep"] == CACHED_SPECIES
    assert result["items"][0]["match_type"] == "id_resolve"


def test_fetch_eggnog_from_cache(store: IndexStore) -> None:
    eggnog = CACHE_DIR / f"{CACHED_SPECIES}_eggNOG.tsv"
    if not eggnog.is_file():
        pytest.skip(f"本地缓存不存在: {eggnog}")

    result = fetch_annotations(
        store,
        species_rep=CACHED_SPECIES,
        biome=BIOME,
        release=RELEASE,
        roles=["eggnog_tsv"],
        convert_gtf=False,
        preview_rows=3,
    )
    assert result["species_rep"] == CACHED_SPECIES
    assert result["source"] == "cache"
    assert any(f["role"] == "eggnog_tsv" and f["status"] == "ready" for f in result["files"])
    assert "eggnog_tsv" in result.get("preview", {})
    preview = result["preview"]["eggnog_tsv"]
    assert len(preview["rows"]) <= 3
    assert preview["columns"]


@patch("tools.mgnify.fetch.list_files", return_value=[])
def test_fetch_remote_skip_when_no_cache(mock_list, store: IndexStore) -> None:
    eggnog = CACHE_DIR / f"{CACHED_SPECIES}_eggNOG.tsv"
    if eggnog.is_file():
        pytest.skip("已有缓存，跳过远程 mock 用例")

    result = fetch_annotations(
        store,
        species_rep="MGYG000000001",
        biome=BIOME,
        release=RELEASE,
        roles=["eggnog_tsv"],
        convert_gtf=False,
        preview_rows=0,
    )
    assert result["files"][0]["status"] == "missing"
    mock_list.assert_called_once()
