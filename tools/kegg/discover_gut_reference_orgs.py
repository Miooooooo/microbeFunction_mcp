"""
Discover ~N gut-associated prokaryote KEGG org codes from list/organism.

Writes kegg_mcp/data/reference_orgs_gut_v{version}.json
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import requests

BASE_URL = "https://rest.kegg.jp"
MIN_INTERVAL = 0.34
DATA_DIR = Path(__file__).resolve().parent / "data"

GUT_GENERA = (
    "Bacteroides",
    "Prevotella",
    "Segatella",
    "Parabacteroides",
    "Phocaeicola",
    "Faecalibacterium",
    "Faecalibacillus",
    "Bifidobacterium",
    "Lactobacillus",
    "Lactiplantibacillus",
    "Ligilactobacillus",
    "Limosilactobacillus",
    "Lacticaseibacillus",
    "Latilactobacillus",
    "Levilactobacillus",
    "Clostridium",
    "Clostridioides",
    "Roseburia",
    "Eubacterium",
    "Agathobacter",
    "Anaerobutyricum",
    "Akkermansia",
    "Escherichia",
    "Streptococcus",
    "Enterococcus",
    "Ruminococcus",
    "Blautia",
    "Megasphaera",
    "Veillonella",
    "Collinsella",
    "Bilophila",
    "Desulfovibrio",
    "Muribaculum",
    "Paramuribaculum",
    "Alistipes",
    "Coprococcus",
    "Eggerthella",
    "Anaerostipes",
    "Hungatella",
    "Butyrivibrio",
    "Lachnospira",
    "Dialister",
    "Fusobacterium",
    "Barnesiella",
    "Odoribacter",
    "Porphyromonas",
    "Selenomonas",
    "Methanobrevibacter",
    "Methanomassiliicoccus",
    "Christensenella",
    "Oscillibacter",
    "Flavonifractor",
    "Intestinimonas",
    "Terrisporobacter",
    "Acidaminococcus",
    "Megamonas",
    "Sutterella",
    "Dorea",
    "Holdemanella",
    "Gemella",
    "Rothia",
    "Actinomyces",
    "Klebsiella",
    "Enterobacter",
    "Citrobacter",
    "Raoultella",
    "Klebsiella",
    "Haemophilus",
    "Campylobacter",
    "Helicobacter",
    "Succinivibrio",
    "Phascolarctobacterium",
    "Anaerotruncus",
    "Butyricicoccus",
    "Tyzzerella",
    "Massiliimalia",
    "Hafnia",
    "Proteus",
)

NAME_KEYWORDS = (
    "intestinal",
    "intestinale",
    "intestin",
    "colon",
    "fecal",
    "faecal",
    "caecal",
    "cecal",
    "rectale",
    "rectalis",
    "gut ",
    " gut",
)

EXCLUDE_NAME = (
    "orangutan",
    "mouse",
    "human",
    "zebra",
    "chicken",
    "pig ",
    "cattle",
    "fish",
    "plant",
    "yeast",
    "virus",
    "plasmid",
)


def fetch_organism_list() -> list[tuple[str, str, str]]:
    time.sleep(MIN_INTERVAL)
    text = requests.get(f"{BASE_URL}/list/organism", timeout=180).text
    rows: list[tuple[str, str, str]] = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        code = parts[1].strip()
        name = parts[2].strip()
        tax = parts[3] if len(parts) > 3 else ""
        if len(code) != 3:
            continue
        if "Bacteria" not in tax and "Archaea" not in tax:
            continue
        rows.append((code, name, tax))
    return rows


def score_candidate(name: str, tax: str) -> int:
    lower = name.lower()
    if any(x in lower for x in EXCLUDE_NAME):
        return -1
    score = 0
    for g in GUT_GENERA:
        if name.startswith(g + " "):
            score += 10
            break
    for kw in NAME_KEYWORDS:
        if kw in lower:
            score += 5
    if "intestinal microbiota" in tax.lower():
        score += 3
    return score


def load_prior_codes(path: Path) -> list[tuple[str, str]]:
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    return [(o["code"], o["name"]) for o in data["organisms"]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=int, default=100)
    parser.add_argument("--version", default="0.1.2")
    parser.add_argument(
        "--seed",
        type=Path,
        default=DATA_DIR / "reference_orgs_gut_v0.1.1.json",
    )
    args = parser.parse_args()

    prior = load_prior_codes(args.seed)
    prior_codes = {c for c, _ in prior}
    rows = fetch_organism_list()

    by_code: dict[str, tuple[str, str]] = {c: (c, n) for c, n in prior}

    ranked: list[tuple[int, str, str]] = []
    for code, name, tax in rows:
        if code in prior_codes:
            continue
        sc = score_candidate(name, tax)
        if sc > 0:
            ranked.append((sc, code, name))
    ranked.sort(key=lambda x: (-x[0], x[1]))

    for _sc, code, name in ranked:
        if len(by_code) >= args.target:
            break
        by_code[code] = (code, name)

  # fill if still short: any Bacteroidota/Firmicutes genus from metadata
    if len(by_code) < args.target:
        for code, name, tax in rows:
            if len(by_code) >= args.target:
                break
            if code in by_code:
                continue
            if any(
                p in tax
                for p in (
                    "Bacteroidota",
                    "Bacillota",
                    "Firmicutes",
                    "Actinomycetota",
                    "Verrucomicrobiota",
                    "Proteobacteria",
                )
            ) and score_candidate(name, tax) >= 0:
                if name.split()[0] in GUT_GENERA or "intestinal" in name.lower():
                    by_code[code] = (code, name)

    organisms = [{"code": c, "name": n} for c, n in sorted(by_code.values(), key=lambda x: x[0])]
    out = {
        "version": args.version,
        "description": f"{len(organisms)} gut-associated prokaryote KEGG codes (auto-discovered + v0.1.1 seed).",
        "target_count": args.target,
        "seed_version": args.seed.name,
        "organisms": organisms,
    }
    out_path = DATA_DIR / f"reference_orgs_gut_v{args.version}.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {len(organisms)} organisms -> {out_path}")


if __name__ == "__main__":
    main()
