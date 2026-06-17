"""
爬取 EBI MGnify Genomes FTP 目录下的子文件夹名称。

目标 URL:
https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/
"""

from __future__ import annotations

import re
import urllib.error
import urllib.request

DEFAULT_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/"
)

# Apache 目录列表中，文件夹行带有 folder.gif 图标
_FOLDER_LINK_RE = re.compile(
    r'<img src="/icons/folder\.gif"[^>]*></td><td><a href="([^"]+)/">',
    re.IGNORECASE,
)

# Apache 目录列表中的文件行（text / unknown 等图标）
_FILE_LINK_RE = re.compile(
    r'<img src="/icons/(?:unknown|text|compressed)\.gif"[^>]*></td><td><a href="([^"]+)">',
    re.IGNORECASE,
)

_USER_AGENT = "Mozilla/5.0 (compatible; mgnify-folder-lister/1.0)"


def _fetch_directory_html(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode(
            resp.headers.get_content_charset() or "utf-8", errors="replace"
        )


def list_subfolders(url: str, timeout: int = 60) -> list[str]:
    """获取任意目录索引页面下的子文件夹名称（按字母序）。"""
    html = _fetch_directory_html(url, timeout=timeout)
    return sorted(_FOLDER_LINK_RE.findall(html))


def list_files(url: str, timeout: int = 60) -> list[str]:
    """获取任意目录索引页面下的文件名（按字母序）。"""
    html = _fetch_directory_html(url, timeout=timeout)
    return sorted(_FILE_LINK_RE.findall(html))


def list_folders(url: str = DEFAULT_URL, timeout: int = 60) -> list[str]:
    """
    获取指定 URL 下的子文件夹名称列表（不含尾部斜杠）。

    Args:
        url: 目录索引页面的 HTTPS 地址。
        timeout: 请求超时秒数。

    Returns:
        按字母序排序的文件夹名称列表。

    Raises:
        urllib.error.URLError: 网络或 HTTP 请求失败。
    """
    return list_subfolders(url, timeout=timeout)


if __name__ == "__main__":
    try:
        names = list_folders()
    except urllib.error.URLError as e:
        raise SystemExit(f"请求失败: {e}") from e

    print(f"共 {len(names)} 个文件夹:\n")
    for name in names:
        print(name)

    print("\nPython list:")
    print(names)
