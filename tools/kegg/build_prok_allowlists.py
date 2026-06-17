"""
Build prokaryote-focused KEGG module / pathway allowlists from KEGG REST.

No gut prevalence — union of terms linked in reference prokaryote genomes,
plus keyword / ID-pattern exclusions for obvious eukaryote-only content.

Usage:
  # Legacy: single prokaryote allowlist from gut reference panel
  python kegg_mcp/build_prok_allowlists.py --version 0.1
  python kegg_mcp/build_prok_allowlists.py --version 0.1.1
  python kegg_mcp/build_prok_allowlists.py --version 0.1.2
  python kegg_mcp/build_prok_allowlists.py --compare

  # Split-domain: separate Bacteria/Archaea MODULE + PATHWAY allowlists
  # built from br08601 taxonomy via stratified per-phylum panel
  python kegg_mcp/build_prok_allowlists.py --split-domains --version 2.0
  python kegg_mcp/build_prok_allowlists.py --split-domains --version 2.0 --no-cache-taxonomy
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

# v0.1 panel (12). Note: mrs is NOT Muribaculum; use pary in v0.1.1+.
REFERENCE_ORGS_V01: list[tuple[str, str]] = [
    ("eco", "Escherichia coli K-12 MG1655"),
    ("bth", "Bacteroides thetaiotaomicron"),
    ("bfr", "Bacteroides fragilis"),
    ("ppu", "Prevotella-class genome (KEGG ppu)"),
    ("pdi", "Parabacteroides distasonis"),
    ("fpr", "Faecalibacterium prausnitzii"),
    ("bli", "Bifidobacterium longum"),
    ("lpl", "Lactiplantibacillus plantarum"),
    ("cpe", "Clostridium perfringens"),
    ("bsu", "Bacillus subtilis"),
    ("lmo", "Listeria monocytogenes"),
    ("mrs", "KEGG mrs (not Muribaculum; kept for v0.1 reproducibility)"),
]

PATHWAY_EXCLUDE_IDS = frozenset({"ko01100", "ko01110", "ko01120"})

PATHWAY_EXCLUDE_PREFIXES = (
    "ko041",
    "ko042",
    "ko043",
    "ko045",
    "ko046",
    "ko047",
    "ko048",
    "ko049",
    "ko05",
    "ko053",
    "ko054",
)

MODULE_EXCLUDE_KEYWORDS = (
    "spliceosome",
    "proteasome",
    "lysosome",
    "cell cycle",
    "mitosis",
    "meiosis",
    "steroid hormone",
    "insulin",
    "glucagon signaling",
    "RNA transport, eukaryotic",
    "nucleotide excision repair, eukaryotic",
    "homologous recombination, eukaryotic",
    "mismatch repair, eukaryotic",
    "Fanconi anemia",
    "p53 signaling",
    "Wnt signaling",
    "Notch signaling",
    "TGF-beta signaling",
    "VEGF signaling",
    "ErbB signaling",
    "JAK-STAT signaling",
    "chemokine signaling",
    "complement and coagulation",
    "Alzheimer disease",
    "Parkinson disease",
    "Huntington disease",
    "amyotrophic lateral sclerosis",
    "pathways of neurodegeneration",
    "human ",
    "coronavirus disease",
    "influenza A",
    "hepatitis ",
    "measles",
    "malaria",
)

PATHWAY_EXCLUDE_KEYWORDS = (
    "human ",
    "Alzheimer",
    "Parkinson",
    "Huntington",
    "amyotrophic lateral sclerosis",
    "pathways of neurodegeneration",
    "coronavirus",
    "influenza A",
    "hepatitis C",
    "measles",
    "malaria",
    "cancer:",
    "colorectal cancer",
    "breast cancer",
    "prostate cancer",
    "melanoma",
    "acute myeloid leukemia",
    "chronic myeloid leukemia",
    "renal cell carcinoma",
    "pancreatic cancer",
    "gastric cancer",
    "hepatocellular carcinoma",
    "non-small cell lung cancer",
    "small cell lung cancer",
    "thyroid cancer",
    "endometrial cancer",
    "basal cell carcinoma",
    "viral carcinogenesis",
    "chemical carcinogenesis",
    "spliceosome",
    "proteasome",
    "lysosome",
    "peroxisome",
)

MODULE_ID_RE = re.compile(r"M\d{5}")

# BRITE organism taxonomy tree (br08601) constants
BRITE_ORGANISMS = "br08601"
TAXONOMY_CACHE = DATA_DIR / "kegg_organisms_taxonomy.json"
BRITE_LEVEL_RE = re.compile(r"^([A-E])(.*)$")
ORG_CODE_RE = re.compile(r"^([a-z]{3,4})\s{2,}(.+)$")


class KeggClient:
    def __init__(self) -> None:
        self._last = 0.0
        self.session = requests.Session()

    def get(self, path: str, retries: int = 4) -> str:
        """GET with throttling + retry on transient network errors.

        400 (bad request, e.g. org with no KEGG data) is raised immediately;
        ConnectionError/Timeout are retried with exponential backoff.
        """
        last_err: Exception | None = None
        for attempt in range(retries):
            elapsed = time.monotonic() - self._last
            if elapsed < MIN_INTERVAL:
                time.sleep(MIN_INTERVAL - elapsed)
            url = f"{BASE_URL}{path}"
            try:
                resp = self.session.get(url, timeout=120)
                self._last = time.monotonic()
                resp.raise_for_status()
                return resp.text
            except requests.exceptions.HTTPError:
                # 4xx/5xx from server; do not retry, propagate
                self._last = time.monotonic()
                raise
            except (
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
            ) as e:
                # transient network issue -> retry with backoff
                self._last = time.monotonic()
                last_err = e
                wait = 2 ** attempt  # 1, 2, 4, 8 s
                time.sleep(wait)
                continue
        assert last_err is not None
        raise last_err


def load_reference_orgs(version: str) -> list[tuple[str, str]]:
    if version == "0.1":
        return list(REFERENCE_ORGS_V01)
    panel_path = DATA_DIR / f"reference_orgs_gut_v{version}.json"
    if panel_path.exists():
        payload = json.loads(panel_path.read_text(encoding="utf-8"))
        return [(o["code"], o["name"]) for o in payload["organisms"]]
    raise ValueError(
        f"Unknown version {version!r}: add {panel_path.name} or use 0.1"
    )


def parse_module_links(text: str) -> set[str]:
    return set(MODULE_ID_RE.findall(text))


def org_path_to_ko(path_id: str) -> str | None:
    if len(path_id) < 8:
        return None
    digits = path_id[3:8]
    if not digits.isdigit():
        return None
    return f"ko{digits}"


def should_exclude_module(module_id: str, name: str) -> bool:
    lower = name.lower()
    return any(kw.lower() in lower for kw in MODULE_EXCLUDE_KEYWORDS)


def should_exclude_pathway(ko_id: str, name: str) -> bool:
    if ko_id in PATHWAY_EXCLUDE_IDS:
        return True
    if any(ko_id.startswith(p) for p in PATHWAY_EXCLUDE_PREFIXES):
        return True
    lower = name.lower()
    return any(kw.lower() in lower for kw in PATHWAY_EXCLUDE_KEYWORDS)


def fetch_all_module_names(client: KeggClient) -> dict[str, str]:
    names: dict[str, str] = {}
    for line in client.get("/list/module").splitlines():
        if not line.strip():
            continue
        mid, _, name = line.partition("\t")
        names[mid.strip()] = name.strip()
    return names


def build_allowlists(
    client: KeggClient,
    orgs: list[tuple[str, str]],
) -> tuple[dict[str, dict], dict[str, dict]]:
    all_module_names = fetch_all_module_names(client)

    module_sources: dict[str, set[str]] = {}
    pathway_sources: dict[str, set[str]] = {}
    pathway_names: dict[str, str] = {}

    for org, _label in orgs:
        mod_text = client.get(f"/link/module/{org}")
        for mid in parse_module_links(mod_text):
            module_sources.setdefault(mid, set()).add(org)

        for line in client.get(f"/list/pathway/{org}").splitlines():
            if not line.strip():
                continue
            pid, _, pname = line.partition("\t")
            pid = pid.strip()
            ko = org_path_to_ko(pid)
            if not ko:
                continue
            pathway_names.setdefault(ko, pname.strip())
            pathway_sources.setdefault(ko, set()).add(org)

    modules_out: dict[str, dict] = {}
    for mid in sorted(module_sources):
        name = all_module_names.get(mid, "")
        if should_exclude_module(mid, name):
            continue
        modules_out[mid] = {
            "module_id": mid,
            "module_name": name,
            "tier": "A",
            "source_orgs": sorted(module_sources[mid]),
            "source_org_count": len(module_sources[mid]),
        }

    pathways_out: dict[str, dict] = {}
    for ko in sorted(pathway_sources):
        name = pathway_names.get(ko, "")
        if should_exclude_pathway(ko, name):
            continue
        pathways_out[ko] = {
            "pathway_id": ko,
            "pathway_name": name,
            "source_orgs": sorted(pathway_sources[ko]),
            "source_org_count": len(pathway_sources[ko]),
        }

    return modules_out, pathways_out


def fetch_taxonomy(client: KeggClient, use_cache: bool = True) -> dict:
    """Fetch and parse br08601 BRITE organism tree; cache to JSON.

    Returns dict with organisms list (code/name/domain/phylum/class) and
    per-domain phylum counts.
    """
    if use_cache and TAXONOMY_CACHE.exists():
        return json.loads(TAXONOMY_CACHE.read_text(encoding="utf-8"))

    text = client.get(f"/get/br:{BRITE_ORGANISMS}")
    (DATA_DIR / "br08601_raw.keg").write_text(text, encoding="utf-8")

    current_domain = None
    current_phylum = None
    current_class = None
    organisms: list[dict] = []

    for ln in text.splitlines():
        m = BRITE_LEVEL_RE.match(ln)
        if not m:
            continue
        level, rest = m.group(1), m.group(2).rstrip()
        if level == "A":
            # superkingdom row; domain set on B-level
            continue
        if level == "B":
            low = rest.lower()
            if "bacteria" in low:
                current_domain = "Bacteria"
            elif "archaea" in low:
                current_domain = "Archaea"
            else:
                current_domain = rest
            current_phylum = None
            current_class = None
        elif level == "C":
            current_phylum = rest
            current_class = None
        elif level == "D":
            current_class = rest
        elif level == "E":
            om = ORG_CODE_RE.match(rest.strip())
            if om and current_domain in ("Bacteria", "Archaea"):
                organisms.append({
                    "code": om.group(1),
                    "name": om.group(2).strip(),
                    "domain": current_domain,
                    "phylum": current_phylum,
                    "class": current_class,
                })

    bac = [o for o in organisms if o["domain"] == "Bacteria"]
    arc = [o for o in organisms if o["domain"] == "Archaea"]
    from collections import Counter
    out = {
        "source": BRITE_ORGANISMS,
        "bacteria_count": len(bac),
        "archaea_count": len(arc),
        "bacteria_phyla": dict(Counter(o["phylum"] for o in bac)),
        "archaea_phyla": dict(Counter(o["phylum"] for o in arc)),
        "organisms": organisms,
    }
    TAXONOMY_CACHE.write_text(
        json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return out


def build_domain_panels(
    taxonomy: dict,
    bac_big_n: int = 2,
    bac_small_n: int = 1,
    bac_big_threshold: int = 100,
    arc_big_n: int = 2,
    arc_small_n: int = 1,
    arc_big_threshold: int = 20,
) -> tuple[list[dict], list[dict]]:
    """Stratified sample per phylum from taxonomy.

    Large phyla get big_n organisms, small phyla get small_n.
    Returns (bacteria_panel, archaea_panel) lists of {code,name,phylum}.
    """
    from collections import defaultdict

    by_dp: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    for o in taxonomy["organisms"]:
        by_dp[o["domain"]][o["phylum"]].append(o)

    def sample(phyla_map, big_n, small_n, threshold):
        panel: list[dict] = []
        for phylum, orgs in sorted(phyla_map.items(), key=lambda x: -len(x[1])):
            n = big_n if len(orgs) >= threshold else small_n
            for o in orgs[:n]:
                panel.append({
                    "code": o["code"],
                    "name": o["name"],
                    "phylum": o["phylum"],
                })
        return panel

    bac = sample(by_dp["Bacteria"], bac_big_n, bac_small_n, bac_big_threshold)
    arc = sample(by_dp["Archaea"], arc_big_n, arc_small_n, arc_big_threshold)
    return bac, arc


def build_split_allowlists(
    client: KeggClient,
    bac_orgs: list[dict],
    arc_orgs: list[dict],
) -> dict:
    """Build MODULE + PATHWAY allowlists separately for Bacteria and Archaea.

    Returns {
      "bacteria": {"modules": {...}, "pathways": {...}},
      "archaea":  {"modules": {...}, "pathways": {...}},
    }
    """
    all_module_names = fetch_all_module_names(client)

    def build_one(panel: list[dict]) -> tuple[dict, dict]:
        module_sources: dict[str, set[str]] = {}
        pathway_sources: dict[str, set[str]] = {}
        pathway_names: dict[str, str] = {}
        skipped: list[str] = []

        for i, o in enumerate(panel, 1):
            org = o["code"]
            try:
                mod_text = client.get(f"/link/module/{org}")
                path_text = client.get(f"/list/pathway/{org}")
            except requests.exceptions.HTTPError as e:
                if e.response is not None and e.response.status_code == 400:
                    skipped.append(org)
                    continue
                raise
            for mid in parse_module_links(mod_text):
                module_sources.setdefault(mid, set()).add(org)
            for line in path_text.splitlines():
                if not line.strip():
                    continue
                pid, _, pname = line.partition("\t")
                pid = pid.strip()
                ko = org_path_to_ko(pid)
                if not ko:
                    continue
                pathway_names.setdefault(ko, pname.strip())
                pathway_sources.setdefault(ko, set()).add(org)
            if i % 10 == 0:
                print(f"    ...{i}/{len(panel)} orgs processed")

        if skipped:
            print(f"    skipped {len(skipped)} orgs with no KEGG data: {', '.join(skipped)}")

        modules_out: dict[str, dict] = {}
        for mid in sorted(module_sources):
            name = all_module_names.get(mid, "")
            if should_exclude_module(mid, name):
                continue
            modules_out[mid] = {
                "module_id": mid,
                "module_name": name,
                "tier": "A",
                "source_orgs": sorted(module_sources[mid]),
                "source_org_count": len(module_sources[mid]),
            }

        pathways_out: dict[str, dict] = {}
        for ko in sorted(pathway_sources):
            name = pathway_names.get(ko, "")
            if should_exclude_pathway(ko, name):
                continue
            pathways_out[ko] = {
                "pathway_id": ko,
                "pathway_name": name,
                "source_orgs": sorted(pathway_sources[ko]),
                "source_org_count": len(pathway_sources[ko]),
            }
        return modules_out, pathways_out

    bac_mod, bac_path = build_one(bac_orgs)
    arc_mod, arc_path = build_one(arc_orgs)
    return {
        "bacteria": {"modules": bac_mod, "pathways": bac_path},
        "archaea": {"modules": arc_mod, "pathways": arc_path},
    }


def write_split_outputs(
    out_dir: Path,
    version: str,
    split_data: dict,
    bac_panel: list[dict],
    arc_panel: list[dict],
) -> None:
    """Write 4 domain-specific allowlists (Bacteria/Archaea x MODULE/PATHWAY) + meta."""
    out_dir.mkdir(parents=True, exist_ok=True)

    def write_module_tsv(path: Path, modules: dict) -> None:
        with path.open("w", encoding="utf-8", newline="") as f:
            f.write("module_id\tmodule_name\ttier\tsource_org_count\tsource_orgs\n")
            for row in modules.values():
                f.write(
                    f"{row['module_id']}\t{row['module_name']}\t{row['tier']}\t"
                    f"{row['source_org_count']}\t{','.join(row['source_orgs'])}\n"
                )

    def write_pathway_files(prefix: str, ko_label: str, pathways: dict) -> tuple[Path, Path]:
        txt = out_dir / f"{prefix}.txt"
        with txt.open("w", encoding="utf-8") as f:
            f.write(f"# {ko_label} ko pathway allowlist v{version}\n")
            f.write("# One ko##### per line; map##### excluded; ko01100/1110/1120 excluded\n")
            for ko in sorted(pathways):
                f.write(f"{ko}\n")
        tsv = out_dir / f"{prefix}.tsv"
        with tsv.open("w", encoding="utf-8", newline="") as f:
            f.write("pathway_id\tpathway_name\tsource_org_count\tsource_orgs\n")
            for row in pathways.values():
                f.write(
                    f"{row['pathway_id']}\t{row['pathway_name']}\t"
                    f"{row['source_org_count']}\t{','.join(row['source_orgs'])}\n"
                )
        return txt, tsv

    bac = split_data["bacteria"]
    arc = split_data["archaea"]

    write_module_tsv(out_dir / f"kegg_module_allowlist_bacteria_v{version}.tsv", bac["modules"])
    write_module_tsv(out_dir / f"kegg_module_allowlist_archaea_v{version}.tsv", arc["modules"])

    write_pathway_files(
        f"prokaryote_ko_pathway_allowlist_bacteria_v{version}", "bacteria", bac["pathways"]
    )
    write_pathway_files(
        f"prokaryote_ko_pathway_allowlist_archaea_v{version}", "archaea", arc["pathways"]
    )

    meta = {
        "version": version,
        "method": "stratified per-phylum panel + KEGG link/module & list/pathway, split by domain",
        "taxonomy_source": BRITE_ORGANISMS,
        "bacteria_panel_count": len(bac_panel),
        "archaea_panel_count": len(arc_panel),
        "bacteria_panel": bac_panel,
        "archaea_panel": arc_panel,
        "bacteria": {
            "module_count": len(bac["modules"]),
            "pathway_count": len(bac["pathways"]),
        },
        "archaea": {
            "module_count": len(arc["modules"]),
            "pathway_count": len(arc["pathways"]),
        },
        "module_exclude_keywords": list(MODULE_EXCLUDE_KEYWORDS),
        "pathway_exclude_keywords": list(PATHWAY_EXCLUDE_KEYWORDS),
        "pathway_excluded_ids": sorted(PATHWAY_EXCLUDE_IDS),
        "pathway_excluded_prefixes": list(PATHWAY_EXCLUDE_PREFIXES),
        "files": {
            "bacteria_modules": f"kegg_module_allowlist_bacteria_v{version}.tsv",
            "archaea_modules": f"kegg_module_allowlist_archaea_v{version}.tsv",
            "bacteria_pathways_txt": f"prokaryote_ko_pathway_allowlist_bacteria_v{version}.txt",
            "bacteria_pathways_tsv": f"prokaryote_ko_pathway_allowlist_bacteria_v{version}.tsv",
            "archaea_pathways_txt": f"prokaryote_ko_pathway_allowlist_archaea_v{version}.txt",
            "archaea_pathways_tsv": f"prokaryote_ko_pathway_allowlist_archaea_v{version}.tsv",
        },
    }
    meta_path = out_dir / f"prokaryote_allowlist_split_meta_v{version}.json"
    meta_path.write_text(
        json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def write_outputs(
    out_dir: Path,
    version: str,
    modules: dict[str, dict],
    pathways: dict[str, dict],
    orgs: list[tuple[str, str]],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    vtag = version.replace(".", "_")

    module_tsv = out_dir / f"kegg_module_allowlist_prok_v{version}.tsv"
    with module_tsv.open("w", encoding="utf-8", newline="") as f:
        f.write("module_id\tmodule_name\ttier\tsource_org_count\tsource_orgs\n")
        for row in modules.values():
            f.write(
                f"{row['module_id']}\t{row['module_name']}\t{row['tier']}\t"
                f"{row['source_org_count']}\t{','.join(row['source_orgs'])}\n"
            )

    pathway_txt = out_dir / f"prokaryote_ko_pathway_allowlist_v{version}.txt"
    with pathway_txt.open("w", encoding="utf-8") as f:
        f.write(f"# prokaryote ko pathway allowlist v{version}\n")
        f.write("# One ko##### per line; map##### excluded; ko01100/1110/1120 excluded\n")
        for ko in sorted(pathways):
            f.write(f"{ko}\n")

    pathway_tsv = out_dir / f"prokaryote_ko_pathway_allowlist_v{version}.tsv"
    with pathway_tsv.open("w", encoding="utf-8", newline="") as f:
        f.write("pathway_id\tpathway_name\tsource_org_count\tsource_orgs\n")
        for row in pathways.values():
            f.write(
                f"{row['pathway_id']}\t{row['pathway_name']}\t"
                f"{row['source_org_count']}\t{','.join(row['source_orgs'])}\n"
            )

    meta = {
        "version": version,
        "method": "union of KEGG link/module and list/pathway over reference prokaryote orgs",
        "prevalence_filter": False,
        "reference_org_count": len(orgs),
        "reference_orgs": [{"code": c, "name": n} for c, n in orgs],
        "module_count": len(modules),
        "pathway_count": len(pathways),
        "pathway_excluded_ids": sorted(PATHWAY_EXCLUDE_IDS),
        "pathway_excluded_prefixes": list(PATHWAY_EXCLUDE_PREFIXES),
        "module_exclude_keywords": list(MODULE_EXCLUDE_KEYWORDS),
        "pathway_exclude_keywords": list(PATHWAY_EXCLUDE_KEYWORDS),
        "files": {
            "modules": module_tsv.name,
            "pathways_txt": pathway_txt.name,
            "pathways_tsv": pathway_tsv.name,
        },
    }
    meta_path = out_dir / f"prokaryote_allowlist_meta_v{version}.json"
    meta_path.write_text(
        json.dumps(meta, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def load_allowlist_ids(path: Path, id_column: str) -> set[str]:
    if not path.exists():
        return set()
    if path.suffix == ".txt":
        return {
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        }
    import csv

    ids: set[str] = set()
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            val = row.get(id_column, "").strip()
            if val:
                ids.add(val)
    return ids


def _load_meta(out_dir: Path, version: str) -> dict:
    path = out_dir / f"prokaryote_allowlist_meta_v{version}.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _org_count(meta: dict) -> int:
    return meta.get("reference_org_count") or len(meta.get("reference_orgs", []))


def compare_versions(out_dir: Path) -> dict:
    versions = ["0.1", "0.1.1", "0.1.2"]
    mods: dict[str, set[str]] = {}
    paths: dict[str, set[str]] = {}
    stats: dict[str, dict] = {}

    for ver in versions:
        mods[ver] = load_allowlist_ids(
            out_dir / f"kegg_module_allowlist_prok_v{ver}.tsv", "module_id"
        )
        paths[ver] = load_allowlist_ids(
            out_dir / f"prokaryote_ko_pathway_allowlist_v{ver}.txt", "pathway_id"
        )
        meta = _load_meta(out_dir, ver)
        stats[ver] = {
            "reference_orgs": _org_count(meta) if meta else None,
            "modules": len(mods[ver]),
            "pathways": len(paths[ver]),
        }

    report: dict = {"summary": stats}
    for i in range(1, len(versions)):
        prev, curr = versions[i - 1], versions[i]
        report[f"modules_only_in_v{curr}_vs_v{prev}"] = sorted(mods[curr] - mods[prev])
        report[f"modules_only_in_v{prev}_vs_v{curr}"] = sorted(mods[prev] - mods[curr])
        report[f"pathways_only_in_v{curr}_vs_v{prev}"] = sorted(
            paths[curr] - paths[prev]
        )
        report[f"pathways_only_in_v{prev}_vs_v{curr}"] = sorted(
            paths[prev] - paths[curr]
        )
        report[f"module_overlap_v{prev}_v{curr}"] = len(mods[prev] & mods[curr])
        report[f"pathway_overlap_v{prev}_v{curr}"] = len(paths[prev] & paths[curr])

    report_path = out_dir / "allowlist_compare_v0.1_v0.1.1_v0.1.2.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return report


def print_compare_report(report: dict) -> None:
    stats = report["summary"]
    print("=== Allowlist v0.1 / v0.1.1 / v0.1.2 ===")
    for ver in ("0.1", "0.1.1", "0.1.2"):
        s = stats.get(ver, {})
        if not s.get("modules"):
            print(f"v{ver}: (files missing)")
            continue
        orgs = s.get("reference_orgs", "?")
        print(f"v{ver}: orgs={orgs}  modules={s['modules']}  pathways={s['pathways']}")

    def delta(a: str, b: str) -> None:
        sa, sb = stats[a], stats[b]
        if sa.get("modules") and sb.get("modules"):
            print(
                f"  v{a} -> v{b}: modules +{sb['modules'] - sa['modules']}, "
                f"pathways +{sb['pathways'] - sa['pathways']}"
            )

    delta("0.1", "0.1.1")
    delta("0.1.1", "0.1.2")
    delta("0.1", "0.1.2")

    for pair in (("0.1.1", "0.1.2"), ("0.1", "0.1.2")):
        key = f"modules_only_in_v{pair[1]}_vs_v{pair[0]}"
        new_m = report.get(key, [])
        keyp = f"pathways_only_in_v{pair[1]}_vs_v{pair[0]}"
        new_p = report.get(keyp, [])
        print(f"New modules v{pair[1]} vs v{pair[0]} ({len(new_m)}):", end=" ")
        print(", ".join(new_m[:15]), "..." if len(new_m) > 15 else "")
        print(f"New pathways v{pair[1]} vs v{pair[0]} ({len(new_p)}):", end=" ")
        print(", ".join(new_p[:12]), "..." if len(new_p) > 12 else "")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--version",
        default="0.1.1",
        help="Allowlist version to build (0.1, 0.1.1, 0.1.2)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DATA_DIR,
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Compare v0.1 / v0.1.1 / v0.1.2 files (no KEGG fetch)",
    )
    parser.add_argument(
        "--split-domains",
        action="store_true",
        help="Build separate Bacteria/Archaea allowlists from br08601 taxonomy "
        "(stratified per-phylum panel). Outputs MODULE + PATHWAY for each domain.",
    )
    parser.add_argument(
        "--no-cache-taxonomy",
        action="store_true",
        help="Re-fetch br08601 taxonomy tree instead of using cached JSON.",
    )
    args = parser.parse_args()

    if args.compare:
        report = compare_versions(args.out_dir)
        print_compare_report(report)
        print(f"Full diff -> {args.out_dir / 'allowlist_compare_v0.1_v0.1.1_v0.1.2.json'}")
        return

    if args.split_domains:
        client = KeggClient()
        taxonomy = fetch_taxonomy(client, use_cache=not args.no_cache_taxonomy)
        bac_panel, arc_panel = build_domain_panels(taxonomy)
        print(
            f"taxonomy: Bacteria={taxonomy['bacteria_count']} "
            f"Archaea={taxonomy['archaea_count']}"
        )
        print(
            f"panel: Bacteria={len(bac_panel)} Archaea={len(arc_panel)} "
            f"(fetching module+pathway per org, ~"
            f"{(len(bac_panel)+len(arc_panel))*2} requests)"
        )
        split_data = build_split_allowlists(client, bac_panel, arc_panel)
        write_split_outputs(args.out_dir, args.version, split_data, bac_panel, arc_panel)

        bac = split_data["bacteria"]
        arc = split_data["archaea"]
        print(
            f"\n=== Split-domain allowlists v{args.version} ===\n"
            f"Bacteria: modules={len(bac['modules'])}  pathways={len(bac['pathways'])}\n"
            f"Archaea:  modules={len(arc['modules'])}  pathways={len(arc['pathways'])}"
        )
        print(f"\nFiles in {args.out_dir}:")
        for f in sorted(args.out_dir.glob(f"*_v{args.version}.*")):
            print(f"  {f.name}")
        return

    orgs = load_reference_orgs(args.version)
    client = KeggClient()
    modules, pathways = build_allowlists(client, orgs)
    write_outputs(args.out_dir, args.version, modules, pathways, orgs)

    print(f"reference_orgs: {len(orgs)}")
    print(f"modules:  {len(modules)}")
    print(f"pathways: {len(pathways)}")

    if (args.out_dir / "kegg_module_allowlist_prok_v0.1.tsv").exists():
        report = compare_versions(args.out_dir)
        print()
        print_compare_report(report)


if __name__ == "__main__":
    main()
