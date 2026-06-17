"""MGnify species_index 冒烟测试（三个参考物种）。"""

from __future__ import annotations

from pathlib import Path

import pytest

from tools.mgnify.build_mgnify_index import REFERENCE_SPECIES, search_species_index
from tools.mgnify.download_species_annotations import species_genome_dir_url

REPO_ROOT = Path(__file__).resolve().parents[1]
INDEX_PATH = REPO_ROOT / "data" / "mgnify" / "human-gut_v2.0_species_index.jsonl"
BIOME = "human-gut"
RELEASE = "v2.0"


@pytest.fixture(scope="module")
def species_index_path() -> Path:
    if not INDEX_PATH.is_file():
        pytest.skip(f"索引未构建: {INDEX_PATH}")
    return INDEX_PATH


@pytest.mark.parametrize("species_name,expected_rep", REFERENCE_SPECIES)
def test_reference_species_hit(
    species_index_path: Path,
    species_name: str,
    expected_rep: str,
) -> None:
    hits = search_species_index(species_index_path, BIOME, RELEASE, species_name)
    assert len(hits) >= 1, f"未命中: {species_name}"
    hit = hits[0]
    assert hit["species_rep"] == expected_rep
    assert hit["catalogue_url"] == species_genome_dir_url(
        expected_rep, version=RELEASE, biome=BIOME
    )
