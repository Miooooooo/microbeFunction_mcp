"""
下载 MGnify human-gut 最新版本目录下的 genomes-all_metadata.tsv。

目标路径示例:
https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv
"""

from __future__ import annotations

import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

from .list_mgnify_folders import list_subfolders

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/"
HUMAN_GUT_URL = f"{BASE_URL}human-gut/"
METADATA_FILENAME = "genomes-all_metadata.tsv"

_VERSION_RE = re.compile(r"^v(\d+(?:\.\d+)*)$", re.IGNORECASE)


def _version_key(name: str) -> tuple[int, ...]:
    """将 v2.0.2 等版本名转为可比较的元组 (2, 0, 2)。"""
    m = _VERSION_RE.match(name)
    if not m:
        return ()
    return tuple(int(part) for part in m.group(1).split("."))


def get_latest_version_folder(url: str = HUMAN_GUT_URL, timeout: int = 60) -> str:
    """
    获取 human-gut 下最新版本文件夹名（如 v2.0.2）。

    仅考虑形如 vX.Y 或 vX.Y.Z 的版本目录。
    """
    folders = list_subfolders(url, timeout=timeout)
    version_folders = [name for name in folders if _VERSION_RE.match(name)]
    if not version_folders:
        raise ValueError(f"未在 {url} 找到版本文件夹")
    return max(version_folders, key=_version_key)


def download_human_gut_metadata(
    output_dir: str | Path = ".",
    timeout: int = 60,
) -> Path:
    """
    下载 human-gut 最新版本下的 genomes-all_metadata.tsv。

    Args:
        output_dir: 本地保存目录，默认当前目录。
        timeout: 目录列表请求超时秒数（下载使用更长超时）。

    Returns:
        已保存文件的本地路径。
    """
    version = get_latest_version_folder(timeout=timeout)
    file_url = f"{HUMAN_GUT_URL}{version}/{METADATA_FILENAME}"
    output_path = Path(output_dir) / METADATA_FILENAME
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"最新版本: {version}")
    print(f"下载地址: {file_url}")
    print(f"保存到:   {output_path.resolve()}")

    req = urllib.request.Request(
        file_url,
        headers={"User-Agent": "Mozilla/5.0 (compatible; mgnify-downloader/1.0)"},
    )
    # 文件约 110MB，下载超时设长一些
    with urllib.request.urlopen(req, timeout=600) as resp, output_path.open("wb") as f:
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 1024 * 1024  # 1MB
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if total:
                pct = downloaded * 100 / total
                print(f"\r进度: {pct:5.1f}% ({downloaded // 1024 // 1024} / {total // 1024 // 1024} MB)", end="")
        if total:
            print()

    print("下载完成。")
    return output_path


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    try:
        download_human_gut_metadata(output_dir=out)
    except (urllib.error.URLError, ValueError) as e:
        raise SystemExit(f"失败: {e}") from e
