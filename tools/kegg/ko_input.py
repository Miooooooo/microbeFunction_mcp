"""Load and normalize KO identifiers from annotation tables or strings."""

from __future__ import annotations

import csv
import re
from pathlib import Path


def normalize_kos(ko_list: list[str]) -> set[str]:
    kos: set[str] = set()
    for raw in ko_list:
        for match in re.findall(r"K\d{5}", raw):
            kos.add(match)
    return kos


def load_kos_from_annotation_tsv(
    path: str | Path,
    column: str = "KEGG",
    delimiter: str = "\t",
) -> list[str]:
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"Annotation file not found: {file_path}")

    values: list[str] = []
    with open(file_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None or column not in reader.fieldnames:
            available = ", ".join(reader.fieldnames or [])
            raise ValueError(
                f"Column {column!r} not found in {file_path}. "
                f"Available columns: {available or '(none)'}"
            )
        for row in reader:
            cell = (row.get(column) or "").strip()
            if cell and cell != "-":
                values.append(cell)
    return values


def resolve_kos(
    ko_list: list[str] | None = None,
    annotation_file: str | Path | None = None,
    kegg_column: str = "KEGG",
) -> tuple[set[str], dict]:
    from_annotation: list[str] = []
    from_list: list[str] = list(ko_list) if ko_list else []

    if annotation_file:
        from_annotation = load_kos_from_annotation_tsv(
            annotation_file, column=kegg_column
        )

    raw_values = from_annotation + from_list
    if not raw_values:
        if annotation_file and not ko_list:
            raise ValueError(
                "No KO identifiers found: annotation file has no usable "
                f"values in column {kegg_column!r}"
            )
        raise ValueError(
            "At least one of annotation_file or ko_list must be provided "
            "with non-empty KO content"
        )

    kos = normalize_kos(raw_values)
    if not kos:
        raise ValueError("No valid KO identifiers (K#####) found after normalization")

    metadata = {
        "from_annotation_count": len(from_annotation),
        "from_list_count": len(from_list),
        "unique_ko_count": len(kos),
    }
    return kos, metadata
