"""
Saturation analysis: how many gut reference orgs are needed for a stable allowlist?

1. Rank gut-associated KEGG organisms (discover_gut_reference_orgs heuristics).
2. Fetch link/module + list/pathway per org (cached).
3. Cumulative union at panel sizes -> module/pathway counts.
4. Recommend panel size and write final comprehensive allowlists.

Usage:
  python kegg_mcp/saturate_gut_panel.py --max-fetch 200 --write-final
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path

import requests

from .discover_gut_reference_orgs import (
    DATA_DIR,
    fetch_organism_list,
    load_prior_codes,
    score_candidate,
)
from .build_prok_allowlists import (
    MODULE_EXCLUDE_KEYWORDS,
    MODULE_ID_RE,
    PATHWAY_EXCLUDE_IDS,
    PATHWAY_EXCLUDE_KEYWORDS,
    PATHWAY_EXCLUDE_PREFIXES,
    org_path_to_ko,
    should_exclude_module,
    should_exclude_pathway,
)

BASE_URL = "https://rest.kegg.jp"
MIN_INTERVAL = 0.34
CACHE_PATH = DATA_DIR / "gut_org_kegg_cache.json"


class KeggClient:
    def __init__(self) -> None:
        self._last = 0.0
        self.session = requests.Session()

    def get(self, path: str) -> str:
        elapsed = time.monotonic() - self._last
        if elapsed < MIN_INTERVAL:
            time.sleep(MIN_INTERVAL - elapsed)
        url = f"{BASE_URL}{path}"
        resp = self.session.get(url, timeout=120)
        self._last = time.monotonic()
        resp.raise_for_status()
        return resp.text


def rank_gut_organisms(seed_path: Path, dedupe_species: bool) -> list[dict]:
    """Return organisms sorted by gut relevance score (high first)."""
    seed = load_prior_codes(seed_path)
    seed_codes = {c for c, _ in seed}
    seed_entries = [
        {"code": c, "name": n, "score": 100, "from_seed": True}
        for c, n in seed
    ]

    rows = fetch_organism_list()
    seen_species: set[str] = set()
    ranked: list[dict] = []

    for code, name, tax in rows:
        if code in seed_codes:
            continue
        sc = score_candidate(name, tax)
        if sc <= 0:
            continue
        if dedupe_species:
            parts = name.split()
            species_key = " ".join(parts[:2]).lower() if len(parts) >= 2 else name.lower()
            if species_key in seen_species:
                continue
            seen_species.add(species_key)
        ranked.append(
            {"code": code, "name": name, "score": sc, "from_seed": False}
        )

    ranked.sort(key=lambda x: (-x["score"], x["code"]))
    return seed_entries + ranked


def fetch_org_annotations(client: KeggClient, code: str) -> dict:
    modules = set(MODULE_ID_RE.findall(client.get(f"/link/module/{code}")))
    pathways: set[str] = set()
    for line in client.get(f"/list/pathway/{code}").splitlines():
        if not line.strip():
            continue
        pid = line.split("\t", 1)[0].strip()
        ko = org_path_to_ko(pid)
        if ko:
            pathways.add(ko)
    return {"modules": sorted(modules), "pathways": sorted(pathways)}


def load_cache() -> dict:
    if CACHE_PATH.exists():
        return json.loads(CACHE_PATH.read_text(encoding="utf-8"))
    return {"orgs": {}}


def save_cache(cache: dict) -> None:
    CACHE_PATH.write_text(json.dumps(cache, indent=2), encoding="utf-8")


def fetch_up_to(
    client: KeggClient,
    ranked: list[dict],
    max_fetch: int,
    cache: dict,
) -> list[dict]:
    org_cache = cache.setdefault("orgs", {})
    fetched: list[dict] = []

    for entry in ranked:
        if len(fetched) >= max_fetch:
            break
        code = entry["code"]
        if code in org_cache:
            rec = {**entry, **org_cache[code]}
            fetched.append(rec)
            continue
        try:
            ann = fetch_org_annotations(client, code)
            rec = {**entry, **ann, "ok": True}
        except requests.HTTPError:
            rec = {**entry, "modules": [], "pathways": [], "ok": False}
        org_cache[code] = {
            "name": entry["name"],
            "modules": rec["modules"],
            "pathways": rec["pathways"],
            "ok": rec.get("ok", True),
            "score": entry["score"],
        }
        fetched.append(rec)
        if len(fetched) % 25 == 0:
            save_cache(cache)
            print(f"  cached {len(fetched)} orgs ...")

    save_cache(cache)
    return fetched


def load_module_names(client: KeggClient) -> dict[str, str]:
    names: dict[str, str] = {}
    for line in client.get("/list/module").splitlines():
        if not line.strip():
            continue
        mid, _, name = line.partition("\t")
        names[mid.strip()] = name.strip()
    return names


def cumulative_allowlist(
    fetched: list[dict],
    n: int,
    module_names: dict[str, str],
) -> tuple[set[str], set[str]]:
    subset = [o for o in fetched[:n] if o.get("ok", True)]
    mod_src: dict[str, set[str]] = {}
    path_src: dict[str, set[str]] = {}

    for org in subset:
        code = org["code"]
        for mid in org.get("modules", []):
            mod_src.setdefault(mid, set()).add(code)
        for ko in org.get("pathways", []):
            path_src.setdefault(ko, set()).add(code)

    modules: set[str] = set()
    for mid in mod_src:
        if should_exclude_module(mid, module_names.get(mid, "")):
            continue
        modules.add(mid)

    pathways: set[str] = set()
    for ko in path_src:
        if should_exclude_pathway(ko, ""):
            continue
        pathways.add(ko)

    return modules, pathways


def recommend_panel_size(rows: list[dict]) -> int:
    """First N where +step gains <= 1 module and <= 1 pathway vs previous checkpoint."""
    if not rows:
        return 0
    best = rows[0]["n_orgs"]
    for i in range(1, len(rows)):
        prev, curr = rows[i - 1], rows[i]
        d_mod = curr["modules"] - prev["modules"]
        d_path = curr["pathways"] - prev["pathways"]
        if d_mod <= 1 and d_path <= 1:
            return curr["n_orgs"]
        best = curr["n_orgs"]
    return best


def write_final_allowlists(
    fetched: list[dict],
    n: int,
    module_names: dict[str, str],
    version: str,
) -> None:
    subset = fetched[:n]
    mod_src: dict[str, set[str]] = {}
    path_src: dict[str, set[str]] = {}

    for org in subset:
        if not org.get("ok", True):
            continue
        code = org["code"]
        for mid in org.get("modules", []):
            mod_src.setdefault(mid, set()).add(code)
        for ko in org.get("pathways", []):
            path_src.setdefault(ko, set()).add(code)

    mod_rows = []
    for mid in sorted(mod_src):
        name = module_names.get(mid, "")
        if should_exclude_module(mid, name):
            continue
        mod_rows.append(
            {
                "module_id": mid,
                "module_name": name,
                "tier": "A",
                "source_org_count": len(mod_src[mid]),
                "source_orgs": ",".join(sorted(mod_src[mid])),
            }
        )

    path_rows = []
    for ko in sorted(path_src):
        if should_exclude_pathway(ko, ""):
            continue
        path_rows.append(
            {
                "pathway_id": ko,
                "source_org_count": len(path_src[ko]),
                "source_orgs": ",".join(sorted(path_src[ko])),
            }
        )

    mod_path = DATA_DIR / f"kegg_module_allowlist_prok_v{version}.tsv"
    with mod_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "module_id",
                "module_name",
                "tier",
                "source_org_count",
                "source_orgs",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(mod_rows)

    txt_path = DATA_DIR / f"prokaryote_ko_pathway_allowlist_v{version}.txt"
    with txt_path.open("w", encoding="utf-8") as f:
        f.write(f"# comprehensive gut prokaryote pathway allowlist v{version}\n")
        for row in path_rows:
            f.write(row["pathway_id"] + "\n")

    panel_path = DATA_DIR / f"reference_orgs_gut_v{version}.json"
    panel_path.write_text(
        json.dumps(
            {
                "version": version,
                "description": f"Saturation-selected panel ({n} orgs) for comprehensive allowlist",
                "organisms": [
                    {"code": o["code"], "name": o["name"]} for o in subset if o.get("ok", True)
                ],
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )

    meta = {
        "version": version,
        "reference_org_count": len([o for o in subset if o.get("ok", True)]),
        "module_count": len(mod_rows),
        "pathway_count": len(path_rows),
        "saturation_note": "Panel size from saturate_gut_panel.py marginal-gain rule",
    }
    (DATA_DIR / f"prokaryote_allowlist_meta_v{version}.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Final allowlist v{version}: {len(mod_rows)} modules, {len(path_rows)} pathways, {n} orgs")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-fetch", type=int, default=200)
    parser.add_argument(
        "--checkpoints",
        default="12,25,50,75,100,125,150,175,200",
        help="Comma-separated panel sizes to report",
    )
    parser.add_argument("--dedupe-species", action="store_true", default=True)
    parser.add_argument("--write-final", action="store_true")
    parser.add_argument("--final-version", default="1.0")
    parser.add_argument(
        "--seed",
        type=Path,
        default=DATA_DIR / "reference_orgs_gut_v0.1.1.json",
    )
    args = parser.parse_args()

    checkpoints = sorted(
        {int(x) for x in args.checkpoints.split(",") if x.strip()}
    )

    print("Ranking gut organisms ...")
    ranked = rank_gut_organisms(args.seed, dedupe_species=args.dedupe_species)
    print(f"  ranked pool: {len(ranked)} (dedupe_species={args.dedupe_species})")

    cache = load_cache()
    client = KeggClient()
    print(f"Fetching KEGG annotations (max {args.max_fetch}) ...")
    fetched = fetch_up_to(client, ranked, args.max_fetch, cache)
    ok_count = sum(1 for o in fetched if o.get("ok", True))
    print(f"  fetched OK: {ok_count}")

    module_names = load_module_names(client)

    # Global ceiling: keyword-filtered KEGG module definitions
    all_mod_ids = set(module_names)
    global_filtered = {
        m for m in all_mod_ids if not should_exclude_module(m, module_names[m])
    }

    curve: list[dict] = []
    max_n = min(max(checkpoints), len(fetched))
    for n in checkpoints:
        if n > len(fetched):
            continue
        mods, paths = cumulative_allowlist(fetched, n, module_names)
        curve.append(
            {
                "n_orgs": n,
                "modules": len(mods),
                "pathways": len(paths),
                "pct_modules_of_global_filtered": round(
                    100.0 * len(mods) / len(global_filtered), 2
                )
                if global_filtered
                else 0,
            }
        )

    # extend to full fetched if not in checkpoints
    if fetched and (not curve or curve[-1]["n_orgs"] < len(fetched)):
        n = len(fetched)
        mods, paths = cumulative_allowlist(fetched, n, module_names)
        if not curve or curve[-1]["n_orgs"] != n:
            curve.append(
                {
                    "n_orgs": n,
                    "modules": len(mods),
                    "pathways": len(paths),
                    "pct_modules_of_global_filtered": round(
                        100.0 * len(mods) / len(global_filtered), 2
                    )
                    if global_filtered
                    else 0,
                }
            )

    recommended = recommend_panel_size(curve)
    report = {
        "ranked_pool_size": len(ranked),
        "fetched_orgs": len(fetched),
        "global_modules_keyword_filtered": len(global_filtered),
        "recommended_panel_size": recommended,
        "curve": curve,
    }

    report_path = DATA_DIR / "saturation_curve_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    tsv_path = DATA_DIR / "saturation_curve_report.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(curve[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(curve)

    print("\n=== Saturation curve ===")
    print(f"{'N_orgs':>8} {'modules':>8} {'pathways':>8} {'%_global_mod':>12}")
    prev_m = prev_p = 0
    for row in curve:
        dm = row["modules"] - prev_m
        dp = row["pathways"] - prev_p
        prev_m, prev_p = row["modules"], row["pathways"]
        print(
            f"{row['n_orgs']:>8} {row['modules']:>8} {row['pathways']:>8} "
            f"{row['pct_modules_of_global_filtered']:>11.2f}%  (+{dm}/+{dp})"
        )
    print(f"\nGlobal module defs (keyword-filtered): {len(global_filtered)}")
    print(f"Recommended panel size (marginal gain <=1): {recommended}")
    print(f"Reports: {tsv_path}")

    if args.write_final:
        n_final = recommended or (curve[-1]["n_orgs"] if curve else len(fetched))
        write_final_allowlists(fetched, n_final, module_names, args.final_version)


if __name__ == "__main__":
    main()
