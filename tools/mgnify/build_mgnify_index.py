"""
从 genomes-all_metadata.tsv 构建 MGnify release_manifest 与 species_index（JSONL）。

聚合与代表行规则见 schema.md §1.8；物种名解析复用 download_species_annotations。
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from .download_human_gut_metadata import BASE_URL, METADATA_FILENAME
from .download_species_annotations import (
    extract_species_from_lineage,
    mgyg_bucket,
    normalize_species_name,
    species_genome_dir_url,
)

SOURCE_BASE = BASE_URL
DEFAULT_METADATA = Path("genomes-all_metadata.tsv")
DEFAULT_OUTPUT_DIR = Path("data/mgnify")

_FTP_RELEASE_RE = re.compile(
    r"mgnify_genomes/([^/]+)/(v[\d.]+)/",
    re.IGNORECASE,
)

REFERENCE_SPECIES: tuple[tuple[str, str], ...] = (
    ("Enterobacter kobei", "MGYG000000235"),
    ("Enterococcus_D casseliflavus", "MGYG000000238"),
    ("Coprococcus sp900767685", "MGYG000000364"),
)


def extract_genus_from_lineage(lineage: str) -> str | None:
    """从 GTDB 风格 Lineage 中提取属名（紧邻 s__ 前的 g__ 段）。"""
    if ";s__" not in lineage:
        return None
    prefix = lineage.rsplit(";s__", 1)[0]
    if ";g__" not in prefix:
        return None
    return prefix.rsplit(";g__", 1)[-1].strip()


def _parse_float(value: str | None) -> float:
    if value is None or value == "":
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def _row_quality(row: dict[str, str]) -> tuple[float, float]:
    return (
        _parse_float(row.get("Completeness")),
        -_parse_float(row.get("Contamination")),
    )


def _is_better_representative(candidate: dict[str, str], current: dict[str, str], species_rep: str) -> bool:
    cand_genome = candidate.get("Genome", "")
    cur_genome = current.get("Genome", "")
    cand_is_rep = cand_genome == species_rep
    cur_is_rep = cur_genome == species_rep
    if cand_is_rep and not cur_is_rep:
        return True
    if cur_is_rep and not cand_is_rep:
        return False
    return _row_quality(candidate) > _row_quality(current)


def detect_release_from_metadata(
    metadata_path: Path,
    biome: str | None = None,
) -> str:
    """从 TSV 行内 FTP_download 推断 release（取首条匹配行）。"""
    with metadata_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ftp = row.get("FTP_download", "")
            match = _FTP_RELEASE_RE.search(ftp)
            if not match:
                continue
            row_biome, release = match.group(1), match.group(2)
            if biome is None or row_biome == biome:
                return release
    raise ValueError(f"无法从 {metadata_path} 的 FTP_download 推断 release")


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def metadata_fetched_at(metadata_path: Path) -> str:
    mtime = metadata_path.stat().st_mtime
    return datetime.fromtimestamp(mtime, tz=timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def build_release_manifest_entry(
    *,
    biome: str,
    release: str,
    metadata_local_path: Path,
    is_latest: bool = False,
) -> dict:
    release_url = f"{SOURCE_BASE}{biome}/{release}/"
    metadata_file_url = f"{release_url}{METADATA_FILENAME}"
    return {
        "manifest_id": f"{biome}:{release}",
        "biome": biome,
        "release": release,
        "source_base": SOURCE_BASE,
        "release_url": release_url,
        "metadata_file_url": metadata_file_url,
        "metadata_local_path": str(metadata_local_path).replace("\\", "/"),
        "metadata_fetched_at": metadata_fetched_at(metadata_local_path),
        "is_latest": is_latest,
        "available_artifacts": ["genomes_all_metadata"],
    }


@dataclass
class SpeciesBucket:
    genome_count: int = 0
    rep_row: dict[str, str] | None = None
    has_genome_equals_rep: bool = False

    def add_row(self, row: dict[str, str], species_rep: str) -> None:
        self.genome_count += 1
        genome = row.get("Genome", "")
        if genome == species_rep:
            if not self.has_genome_equals_rep:
                self.rep_row = row
                self.has_genome_equals_rep = True
            elif self.rep_row is not None and _is_better_representative(row, self.rep_row, species_rep):
                self.rep_row = row
        elif not self.has_genome_equals_rep:
            if self.rep_row is None:
                self.rep_row = row
            elif _is_better_representative(row, self.rep_row, species_rep):
                self.rep_row = row


def row_to_species_index(
    *,
    biome: str,
    release: str,
    species_rep: str,
    bucket: SpeciesBucket,
    built_at: str,
) -> dict:
    row = bucket.rep_row
    if row is None:
        raise ValueError(f"物种代表 {species_rep} 无代表行")

    lineage = row.get("Lineage", "")
    species_name = extract_species_from_lineage(lineage) or ""
    bucket_id = mgyg_bucket(species_rep)

    return {
        "index_id": f"{biome}:{release}:{species_rep}",
        "species_rep": species_rep,
        "biome": biome,
        "release": release,
        "mgyg_bucket": bucket_id,
        "species_name": species_name,
        "species_name_norm": normalize_species_name(species_name) if species_name else "",
        "genus_name": extract_genus_from_lineage(lineage) or "",
        "lineage": lineage,
        "genome_type": row.get("Genome_type", ""),
        "length_bp": int(_parse_float(row.get("Length"))),
        "completeness": _parse_float(row.get("Completeness")),
        "contamination": _parse_float(row.get("Contamination")),
        "genome_count": bucket.genome_count,
        "catalogue_url": species_genome_dir_url(species_rep, version=release, biome=biome),
        "cache_dir": f"downloads/{species_rep}/genome/",
        "status": "missing",
        "artifacts": {},
        "last_seen_at": built_at,
        "last_checked_at": built_at,
        "last_error": None,
    }


def build_species_index(
    metadata_path: Path,
    biome: str,
    release: str,
) -> list[dict]:
    buckets: dict[str, SpeciesBucket] = {}
    built_at = utc_now_iso()

    with metadata_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species_rep = row.get("Species_rep", "").strip()
            if not species_rep:
                continue
            bucket = buckets.setdefault(species_rep, SpeciesBucket())
            bucket.add_row(row, species_rep)

    entries = [
        row_to_species_index(
            biome=biome,
            release=release,
            species_rep=species_rep,
            bucket=bucket,
            built_at=built_at,
        )
        for species_rep, bucket in sorted(buckets.items())
    ]
    return entries


def merge_manifests(existing: list[dict], new_entry: dict) -> list[dict]:
    by_id = {item["manifest_id"]: item for item in existing if "manifest_id" in item}
    by_id[new_entry["manifest_id"]] = new_entry
    return sorted(by_id.values(), key=lambda x: (x.get("biome", ""), x.get("release", "")))


def write_release_manifest(output_dir: Path, entry: dict) -> Path:
    path = output_dir / "release_manifest.json"
    existing: list[dict] = []
    if path.exists():
        with path.open(encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, list):
            existing = data
        elif isinstance(data, dict):
            existing = [data]
    merged = merge_manifests(existing, entry)
    output_dir.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(merged, f, ensure_ascii=False, indent=2)
        f.write("\n")
    return path


def write_species_index_jsonl(output_dir: Path, biome: str, release: str, entries: list[dict]) -> Path:
    path = output_dir / f"{biome}_{release}_species_index.jsonl"
    output_dir.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for entry in entries:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    return path


def search_species_index(
    jsonl_path: Path,
    biome: str,
    release: str,
    query: str,
) -> list[dict]:
    """
    在 species_index JSONL 中按规范化物种名精确匹配（MVP）。
    """
    target = normalize_species_name(query)
    hits: list[dict] = []
    with jsonl_path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            if entry.get("biome") != biome or entry.get("release") != release:
                continue
            if entry.get("species_name_norm") == target:
                hits.append(entry)
    return hits


def build_mgnify_index(
    *,
    biome: str,
    release: str,
    metadata_path: Path,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    is_latest: bool = False,
) -> tuple[Path, Path, int]:
    manifest_entry = build_release_manifest_entry(
        biome=biome,
        release=release,
        metadata_local_path=metadata_path.resolve(),
        is_latest=is_latest,
    )
    entries = build_species_index(metadata_path, biome, release)
    manifest_path = write_release_manifest(output_dir, manifest_entry)
    index_path = write_species_index_jsonl(output_dir, biome, release, entries)
    return manifest_path, index_path, len(entries)


def _print_reference_hits(index_path: Path, biome: str, release: str) -> None:
    print("\n参考物种命中:")
    for species_name, expected_rep in REFERENCE_SPECIES:
        hits = search_species_index(index_path, biome, release, species_name)
        if not hits:
            print(f"  {species_name}: 未命中")
            continue
        hit = hits[0]
        rep = hit.get("species_rep")
        url = hit.get("catalogue_url")
        expected_url = species_genome_dir_url(expected_rep, version=release, biome=biome)
        ok = rep == expected_rep and url == expected_url
        status = "OK" if ok else "MISMATCH"
        print(f"  [{status}] {species_name}")
        print(f"         species_rep={rep} (期望 {expected_rep})")
        print(f"         catalogue_url={url}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="构建 MGnify release_manifest 与 species_index")
    parser.add_argument("--biome", default="human-gut")
    parser.add_argument("--release", default="v2.0")
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--detect-release-from-metadata",
        action="store_true",
        help="从 TSV 行内 FTP_download 推断 release",
    )
    args = parser.parse_args(argv)

    metadata_path = args.metadata
    if not metadata_path.is_file():
        print(f"metadata 文件不存在: {metadata_path}", file=sys.stderr)
        return 1

    release = args.release
    if args.detect_release_from_metadata:
        release = detect_release_from_metadata(metadata_path, biome=args.biome)
        print(f"从 metadata 推断 release: {release}")

    print(f"构建索引: biome={args.biome} release={release}")
    print(f"metadata: {metadata_path.resolve()}")

    manifest_path, index_path, count = build_mgnify_index(
        biome=args.biome,
        release=release,
        metadata_path=metadata_path,
        output_dir=args.output_dir,
    )

    print(f"release_manifest -> {manifest_path.resolve()}")
    print(f"species_index    -> {index_path.resolve()} ({count} 条)")
    _print_reference_hits(index_path, args.biome, release)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
