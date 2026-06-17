"""
根据菌种名，从 MGnify human-gut 最新版本下载物种代表基因组的注释文件（GFF + 全部 TSV）。

数据存储规律（human-gut/vX.Y.Z）:
  {biome}/{version}/species_catalogue/{MGYG前10位}/{Species_rep}/genome/
    - {Species_rep}.gff
    - {Species_rep}_*.tsv  (InterProScan, eggNOG, KEGG 等)

genomes-all_metadata.tsv 中:
  - Lineage 末尾为物种名，格式 s__{属}_{种} 或 s__{属} {种}
  - Species_rep 列为物种代表基因组 ID（用于 species_catalogue 路径）
  - FTP_download 中 all_genomes/{ID去掉末2位}/{Genome_ID}/ 与上述分桶规则一致
"""

from __future__ import annotations

import csv
import gzip
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .download_human_gut_metadata import (
    BASE_URL,
    HUMAN_GUT_URL,
    get_latest_version_folder,
)
from .list_mgnify_folders import list_files

METADATA_DEFAULT = Path("genomes-all_metadata.tsv")
_USER_AGENT = "Mozilla/5.0 (compatible; mgnify-species-downloader/1.0)"


@dataclass(frozen=True)
class SpeciesHit:
    species_name: str
    species_rep: str
    lineage_example: str
    genome_count: int


def mgyg_bucket(mgyg_id: str) -> str:
    """
    MGnify ID 分桶：去掉末尾 2 位数字。

    例: MGYG000000238 -> MGYG0000002
    与 FTP 路径 all_genomes/{bucket}/{genome_id}/ 一致。
    """
    if not mgyg_id.startswith("MGYG") or len(mgyg_id) < 12:
        raise ValueError(f"无效的 MGnify ID: {mgyg_id}")
    return mgyg_id[:-2]


def extract_species_from_lineage(lineage: str) -> str | None:
    """从 GTDB 风格 Lineage 中提取物种名（s__ 之后）。"""
    if ";s__" not in lineage:
        return None
    return lineage.rsplit(";s__", 1)[-1].strip()


def normalize_species_name(name: str) -> str:
    return " ".join(name.strip().split()).casefold()


def lookup_species_in_metadata(
    species_name: str,
    metadata_path: Path = METADATA_DEFAULT,
) -> SpeciesHit:
    """
    在 genomes-all_metadata.tsv 中按物种名查找 Species_rep。

    匹配 Lineage 中的 s__ 物种名（大小写不敏感，忽略多余空格）。
    """
    target = normalize_species_name(species_name)
    species_rep: str | None = None
    lineage_example = ""
    count = 0

    with metadata_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            lineage = row.get("Lineage", "")
            sp = extract_species_from_lineage(lineage)
            if sp is None or normalize_species_name(sp) != target:
                continue
            count += 1
            rep = row["Species_rep"]
            if species_rep is None:
                species_rep = rep
                lineage_example = lineage
            elif species_rep != rep:
                raise ValueError(
                    f"物种 '{species_name}' 对应多个 Species_rep: "
                    f"{species_rep} 与 {rep}"
                )

    if species_rep is None:
        raise ValueError(f"在 metadata 中未找到物种: {species_name}")

    return SpeciesHit(
        species_name=species_name,
        species_rep=species_rep,
        lineage_example=lineage_example,
        genome_count=count,
    )


def species_genome_dir_url(
    species_rep: str,
    version: str | None = None,
    biome: str = "human-gut",
) -> str:
    version = version or get_latest_version_folder()
    bucket = mgyg_bucket(species_rep)
    return (
        f"{BASE_URL}{biome}/{version}/species_catalogue/"
        f"{bucket}/{species_rep}/genome/"
    )


def annotation_filenames(remote_files: list[str]) -> list[str]:
    """筛选 GFF 与 TSV 注释文件。"""
    selected = []
    for name in remote_files:
        lower = name.lower()
        if lower.endswith(".gff") or lower.endswith(".gff.gz") or lower.endswith(".tsv"):
            selected.append(name)
    return sorted(selected)


def _parse_gff3_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def _feature_id(attrs: dict[str, str], fallback: str) -> str:
    return attrs.get("ID") or attrs.get("locus_tag") or attrs.get("Name") or fallback


def _format_gtf_attributes(gene_id: str, extra: dict[str, str]) -> str:
    """GTF 第 9 列：gene_id / transcript_id 为必需，其余属性保留。"""
    skip = {"ID", "Parent"}
    parts = [
        f'gene_id "{gene_id}"',
        f'transcript_id "{gene_id}"',
    ]
    for key, value in extra.items():
        if key in skip:
            continue
        escaped = value.replace('"', '\\"')
        parts.append(f'{key} "{escaped}"')
    return "; ".join(parts) + ";"


def gff_to_gtf(gff_path: Path, gtf_path: Path | None = None) -> Path:
    """
    将 GFF3 转为 GTF。

    MGnify 原核注释多为 CDS / ncRNA 等扁平特征；为每条记录生成
    gene、transcript、exon（CDS 另加 CDS 行），便于下游工具使用。
    """
    gff_path = Path(gff_path)
    gtf_path = Path(gtf_path or gff_path.with_suffix(".gtf"))

    open_fn = gzip.open if gff_path.suffix == ".gz" else open
    mode = "rt" if gff_path.suffix == ".gz" else "r"

    with open_fn(gff_path, mode, encoding="utf-8") as fin, gtf_path.open(
        "w", encoding="utf-8"
    ) as fout:
        fout.write(f"##gtf-version 2.2\n")
        fout.write(f"##source converted from {gff_path.name}\n")

        for raw in fin:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                if line.startswith("##sequence-region"):
                    fout.write(line + "\n")
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                continue

            seqid, source, feature, start, end, score, strand, frame, attr_str = cols
            attrs = _parse_gff3_attributes(attr_str)
            gene_id = _feature_id(attrs, f"{seqid}_{start}_{end}")
            gtf_attr = _format_gtf_attributes(gene_id, attrs)
            base = [seqid, source, "", start, end, score, strand, frame, gtf_attr]

            for gtf_feature in ("gene", "transcript", "exon"):
                row = base.copy()
                row[2] = gtf_feature
                row[7] = "."
                fout.write("\t".join(row) + "\n")

            if feature == "CDS":
                row = base.copy()
                row[2] = "CDS"
                row[7] = frame
                fout.write("\t".join(row) + "\n")
            else:
                row = base.copy()
                row[2] = feature
                row[7] = "."
                fout.write("\t".join(row) + "\n")

    return gtf_path


def convert_downloaded_gff_to_gtf(paths: list[Path]) -> list[Path]:
    """将路径列表中的 GFF / GFF.GZ 转为同目录下的 .gtf 文件。"""
    gtf_paths: list[Path] = []
    for path in paths:
        lower = path.name.lower()
        if not (lower.endswith(".gff") or lower.endswith(".gff.gz")):
            continue
        if lower.endswith(".gff.gz"):
            gtf_out = path.with_name(path.name.replace(".gff.gz", ".gtf"))
        else:
            gtf_out = path.with_suffix(".gtf")
        print(f"转换 GFF -> GTF: {path.name} -> {gtf_out.name}")
        gtf_paths.append(gff_to_gtf(path, gtf_out))
    return gtf_paths


def _download_file(url: str, dest: Path, timeout: int = 600) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp, dest.open("wb") as f:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def download_species_annotations(
    species_name: str,
    output_dir: str | Path = ".",
    metadata_path: Path = METADATA_DEFAULT,
    version: str | None = None,
    biome: str = "human-gut",
    timeout: int = 60,
    convert_to_gtf: bool = True,
) -> list[Path]:
    """
    按菌种名下载物种代表基因组 genome/ 目录下全部 GFF 与 TSV 文件。

    Args:
        convert_to_gtf: 下载完成后是否将 GFF 转为 GTF（默认为是）。

    Returns:
        已下载/生成的本地文件路径列表（含 GTF）。
    """
    hit = lookup_species_in_metadata(species_name, metadata_path)
    version = version or get_latest_version_folder(f"{BASE_URL}{biome}/", timeout=timeout)
    genome_url = species_genome_dir_url(hit.species_rep, version=version, biome=biome)

    print(f"物种:         {hit.species_name}")
    print(f"Species_rep:  {hit.species_rep}")
    print(f"分桶目录:     {mgyg_bucket(hit.species_rep)}")
    print(f"基因组条目数: {hit.genome_count}")
    print(f"版本:         {version}")
    print(f"远程目录:     {genome_url}")

    remote_files = list_files(genome_url, timeout=timeout)
    to_download = annotation_filenames(remote_files)
    if not to_download:
        raise ValueError(f"目录中未找到 GFF/TSV 文件: {genome_url}")

    out_dir = Path(output_dir) / hit.species_rep / "genome"
    out_dir.mkdir(parents=True, exist_ok=True)

    saved: list[Path] = []
    for i, filename in enumerate(to_download, 1):
        file_url = genome_url + filename
        dest = out_dir / filename
        print(f"[{i}/{len(to_download)}] 下载 {filename} ...")
        _download_file(file_url, dest)
        saved.append(dest)

    if convert_to_gtf:
        gtf_files = convert_downloaded_gff_to_gtf(saved)
        saved.extend(gtf_files)

    print(f"\n完成，共 {len(saved)} 个文件 -> {out_dir.resolve()}")
    return saved


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(
            "用法: python download_species_annotations.py "
            '"Enterococcus_D casseliflavus" [输出目录] [metadata.tsv路径]'
        )

    species = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "./downloads"
    meta = Path(sys.argv[3]) if len(sys.argv) > 3 else METADATA_DEFAULT

    try:
        download_species_annotations(species, output_dir=out, metadata_path=meta)
    except (urllib.error.URLError, ValueError, OSError) as e:
        raise SystemExit(f"失败: {e}") from e
