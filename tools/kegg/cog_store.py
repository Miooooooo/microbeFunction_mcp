"""Local COG catalog loaded from COG.csv for lookup and MCP resources."""

from __future__ import annotations

import csv
import json
import os
import re
from dataclasses import asdict, dataclass
from pathlib import Path

COG_ID_RE = re.compile(r"COG\d{4,5}", re.IGNORECASE)

COG_CATEGORY_NAMES: dict[str, str] = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton",
}


def cog_csv_candidates() -> list[Path]:
    """Search order for the COG catalog file (first existing file wins)."""
    pkg_data = Path(__file__).resolve().parent / "data" / "COG.csv"
    project_root = Path(__file__).resolve().parent.parent / "COG.csv"
    return [pkg_data, project_root]


def resolve_cog_csv_path(explicit: str | Path | None = None) -> Path:
    """
    Resolve COG.csv location.

    Priority: explicit argument > KEGG_MCP_COG_CSV > kegg_mcp/data/COG.csv >
    project-root COG.csv (legacy).
    """
    if explicit is not None:
        path = Path(explicit)
        if not path.is_file():
            raise FileNotFoundError(f"COG catalog not found: {path}")
        return path

    env = os.environ.get("KEGG_MCP_COG_CSV", "").strip()
    if env:
        path = Path(env)
        if not path.is_file():
            raise FileNotFoundError(
                f"KEGG_MCP_COG_CSV points to missing file: {path}"
            )
        return path

    for candidate in cog_csv_candidates():
        if candidate.is_file():
            return candidate

    tried = ", ".join(str(p) for p in cog_csv_candidates())
    raise FileNotFoundError(
        "COG catalog not found. Options: "
        "(1) set env KEGG_MCP_COG_CSV to your COG.csv path; "
        f"(2) place COG.csv at one of: {tried}"
    )


def default_cog_csv_path() -> Path:
    return resolve_cog_csv_path()


def is_cog_catalog_available() -> bool:
    try:
        resolve_cog_csv_path()
        return True
    except FileNotFoundError:
        return False


def normalize_cog_id(value: str) -> str | None:
    if not value or not str(value).strip():
        return None
    text = str(value).strip()
    match = COG_ID_RE.search(text.replace(":", " "))
    if not match:
        return None
    return match.group(0).upper()


def split_cat_codes(cat_code: str) -> list[str]:
    return [c for c in (cat_code or "").upper() if c.isalpha()]


@dataclass(frozen=True)
class CogEntry:
    cog_id: str
    cat_code: str
    annotation: str
    gene_symbol: str
    pathway: str
    pubmed: str
    pdb: str

    def to_dict(self) -> dict:
        cats = split_cat_codes(self.cat_code)
        return {
            **asdict(self),
            "cat_codes": cats,
            "cat_names": [COG_CATEGORY_NAMES.get(c, c) for c in cats],
        }


class CogCatalog:
    """In-memory index over NCBI COG.csv."""

    def __init__(self, csv_path: str | Path | None = None):
        self.csv_path = (
            resolve_cog_csv_path(csv_path) if csv_path is not None else resolve_cog_csv_path()
        )
        self._by_id: dict[str, CogEntry] = {}
        self._loaded = False

    def load(self) -> None:
        if self._loaded:
            return
        with open(self.csv_path, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError(f"COG.csv has no header: {self.csv_path}")
            for row in reader:
                raw_id = (row.get("COG") or row.get("cog") or "").strip()
                cog_id = normalize_cog_id(raw_id)
                if not cog_id:
                    continue
                self._by_id[cog_id] = CogEntry(
                    cog_id=cog_id,
                    cat_code=(row.get("Cat") or row.get("cat") or "").strip(),
                    annotation=(row.get("Annotation") or "").strip(),
                    gene_symbol=(row.get("Gene") or "").strip(),
                    pathway=(row.get("Pathway") or "").strip(),
                    pubmed=(row.get("PubMed") or "").strip(),
                    pdb=(row.get("PDB") or "").strip(),
                )
        if not self._by_id:
            raise ValueError(f"No COG entries loaded from {self.csv_path}")
        self._loaded = True

    def lookup(self, cog_id: str) -> CogEntry | None:
        self.load()
        key = normalize_cog_id(cog_id)
        if not key:
            return None
        return self._by_id.get(key)

    def lookup_many(self, cog_ids: list[str]) -> dict[str, dict | None]:
        self.load()
        out: dict[str, dict | None] = {}
        for raw in cog_ids:
            key = normalize_cog_id(raw)
            if not key:
                continue
            entry = self._by_id.get(key)
            out[key] = entry.to_dict() if entry else None
        return out

    def search(
        self,
        query: str,
        *,
        limit: int = 20,
        cat: str = "",
    ) -> list[dict]:
        self.load()
        q = (query or "").strip().lower()
        if not q:
            return []
        cat_letters = {c.upper() for c in cat if c.isalpha()} if cat else set()
        limit = max(1, min(limit, 200))
        hits: list[dict] = []
        for entry in self._by_id.values():
            if cat_letters and not cat_letters.intersection(split_cat_codes(entry.cat_code)):
                continue
            blob = " ".join(
                [
                    entry.cog_id,
                    entry.cat_code,
                    entry.annotation,
                    entry.gene_symbol,
                    entry.pathway,
                ]
            ).lower()
            if q in blob:
                hits.append(entry.to_dict())
                if len(hits) >= limit:
                    break
        return hits

    def meta(self) -> dict:
        self.load()
        cat_counts: dict[str, int] = {}
        for entry in self._by_id.values():
            for letter in split_cat_codes(entry.cat_code):
                cat_counts[letter] = cat_counts.get(letter, 0) + 1
        return {
            "resource": "cog/catalog",
            "version": 1,
            "source_file": str(self.csv_path),
            "entry_count": len(self._by_id),
            "category_legend": COG_CATEGORY_NAMES,
            "category_entry_counts": dict(sorted(cat_counts.items())),
        }


_catalog: CogCatalog | None = None


def get_cog_catalog(csv_path: str | Path | None = None) -> CogCatalog:
    global _catalog
    if csv_path is not None:
        return CogCatalog(csv_path)
    if _catalog is None:
        _catalog = CogCatalog()
    return _catalog


def entry_as_json(entry: CogEntry | None, cog_id: str) -> str:
    if entry is None:
        return json.dumps(
            {"cog_id": normalize_cog_id(cog_id) or cog_id, "found": False},
            ensure_ascii=False,
            indent=2,
        )
    return json.dumps({"found": True, **entry.to_dict()}, ensure_ascii=False, indent=2)


def meta_as_json() -> str:
    return json.dumps(get_cog_catalog().meta(), ensure_ascii=False, indent=2)
