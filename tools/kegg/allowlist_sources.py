"""
Allowlist / blacklist sources for prokaryote KEGG filtering.

This module centralizes all filter rules (keyword blacklists, ID blacklists,
prefix blacklists) and allowlist metadata as importable, versioned sources.
Both blacklists and allowlists are treated as first-class "sources" that can
be loaded, queried, and written to disk independently of the build pipeline.

Design
------
- ``BlacklistEntry`` dataclass: one row with keyword + category + reason +
  applies_to ("module" | "pathway" | "both").
- ``_BLACKLIST_ENTRIES``: single source of truth, organized by biological
  category for readability.
- Derived constants (``MODULE_EXCLUDE_KEYWORDS``, ``PATHWAY_EXCLUDE_KEYWORDS``,
  ``PATHWAY_EXCLUDE_IDS``, ``PATHWAY_EXCLUDE_PREFIXES``) are computed from the
  table so there is exactly one place to edit.
- ``write_blacklist_files()`` emits TSV + TXT files mirroring the allowlist
  output style (one entry per line, full metadata in TSV).
- ``load_allowlist()`` / ``load_blacklist()`` provide uniform source access.

Usage
-----
    from allowlist_sources import (
        MODULE_EXCLUDE_KEYWORDS,
        PATHWAY_EXCLUDE_KEYWORDS,
        PATHWAY_EXCLUDE_IDS,
        PATHWAY_EXCLUDE_PREFIXES,
        write_blacklist_files,
        load_blacklist,
    )
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"
BLACKLIST_VERSION = "1.0"


# ---------------------------------------------------------------------------
# Blacklist entry schema
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BlacklistEntry:
    """One blacklist rule.

    - keyword: case-insensitive substring matched against KEGG module/pathway
      names. For ID/prefix rules, stored in ``value`` instead via the
      ``rule_type`` field on the wrapping record.
    - category: biological category (organelle, signaling, immune, ...).
    - reason: why this is excluded from prokaryote allowlists.
    - applies_to: "module" | "pathway" | "both".
    """
    keyword: str
    category: str
    reason: str
    applies_to: str


# ---------------------------------------------------------------------------
# Single source of truth: all eukaryote-only / disease / human-specific rules
# ---------------------------------------------------------------------------
# Each tuple: (keyword, category, reason, applies_to)
# applies_to: "module" | "pathway" | "both"
_BLACKLIST_ENTRIES: list[tuple[str, str, str, str]] = [
    # === 1. 细胞器 (Eukaryote-only organelles & transport) ===
    ("spliceosome", "organelle", "真核mRNA剪接体，原核无内含子剪接机制", "both"),
    ("proteasome", "organelle", "真核26S蛋白酶体（原核用Lon/HslUV/Clp系统）", "both"),
    ("lysosome", "organelle", "溶酶体为真核特有细胞器", "both"),
    ("peroxisome", "organelle", "过氧化物酶体为真核特有细胞器", "both"),
    ("mitochondria", "organelle", "线粒体为真核特有细胞器", "both"),
    ("mitochondrial", "organelle", "线粒体相关过程，真核特有", "both"),
    ("chloroplast", "organelle", "叶绿体为真核（植物/藻类）特有", "both"),
    ("endoplasmic reticulum", "organelle", "内质网为真核特有细胞器", "both"),
    ("Golgi", "organelle", "高尔基体为真核特有细胞器", "both"),
    ("nuclear pore", "organelle", "核孔复合体为真核特有", "both"),
    ("nucleocytoplasmic transport", "organelle", "核质运输为真核特有", "both"),
    ("SNARE", "organelle", "SNARE囊泡融合机制为真核特有", "both"),
    ("vesicular transport", "organelle", "囊泡运输为真核特有", "both"),
    ("endocytosis", "organelle", "内吞作用为真核特有", "both"),
    ("exocytosis", "organelle", "胞吐作用为真核特有", "both"),
    ("phagocytosis", "organelle", "吞噬作用为真核特有", "both"),

    # === 2. 细胞骨架与连接 (Cytoskeleton & junctions) ===
    ("actin cytoskeleton", "cytoskeleton_junction", "肌动蛋白骨架调控为真核特有", "both"),
    ("tight junction", "cytoskeleton_junction", "紧密连接为真核特有", "both"),
    ("focal adhesion", "cytoskeleton_junction", "黏着斑为真核特有", "both"),
    ("adherens junction", "cytoskeleton_junction", "黏附连接为真核特有", "both"),
    ("gap junction", "cytoskeleton_junction", "缝隙连接为真核特有", "both"),

    # === 3. 基因表达调控 (Gene expression & DNA repair eukaryote-specific) ===
    ("ribosome biogenesis", "gene_expression", "rRNA加工/核糖体组装为真核特有", "both"),
    ("histone modification", "gene_expression", "组蛋白修饰为真核特有", "both"),
    ("chromatin remodeling", "gene_expression", "染色质重塑为真核特有", "both"),
    ("basal transcription factors", "gene_expression", "真核基础转录因子（TFIID等）", "both"),
    ("RNA transport, eukaryotic", "gene_expression", "真核mRNA核质运输", "module"),
    ("nucleotide excision repair, eukaryotic", "gene_expression", "真核NER模块", "module"),
    ("homologous recombination, eukaryotic", "gene_expression", "真核HR模块", "module"),
    ("mismatch repair, eukaryotic", "gene_expression", "真核MMR模块", "module"),

    # === 4. 细胞周期与分裂 (Cell cycle & division) ===
    ("cell cycle", "cell_cycle", "真核细胞周期（CDK/cyclin）", "both"),
    ("mitosis", "cell_cycle", "有丝分裂为真核特有", "both"),
    ("meiosis", "cell_cycle", "减数分裂为真核特有", "both"),
    ("Fanconi anemia", "cell_cycle", "FA通路为真核DNA损伤应答", "module"),

    # === 5. 信号转导 (Signal transduction) ===
    ("steroid hormone", "signaling", "类固醇激素信号为真核特有", "both"),
    ("insulin", "signaling", "胰岛素信号为动物特有", "both"),
    ("glucagon signaling", "signaling", "胰高血糖素信号为动物特有", "both"),
    ("p53 signaling", "signaling", "p53通路为真核特有", "both"),
    ("Wnt signaling", "signaling", "Wnt通路为动物特有", "both"),
    ("Notch signaling", "signaling", "Notch通路为动物特有", "both"),
    ("TGF-beta signaling", "signaling", "TGF-β通路为动物特有", "both"),
    ("VEGF signaling", "signaling", "VEGF通路为动物特有", "both"),
    ("ErbB signaling", "signaling", "ErbB/EGFR通路为动物特有", "both"),
    ("JAK-STAT signaling", "signaling", "JAK-STAT为动物特有", "both"),
    ("chemokine signaling", "signaling", "趋化因子信号为动物特有", "both"),
    ("Hedgehog signaling", "signaling", "Hedgehog通路为动物特有", "both"),
    ("PI3K-Akt", "signaling", "PI3K-Akt通路为真核特有", "both"),
    ("MAPK signaling", "signaling", "真核MAPK级联（ERK/JNK/p38）", "both"),
    ("AMPK signaling", "signaling", "AMPK为真核能量感受器", "both"),
    ("mTOR signaling", "signaling", "mTOR为真核特有", "both"),
    ("FoxO signaling", "signaling", "FoxO转录因子为动物特有", "both"),
    ("HIF-1 signaling", "signaling", "HIF-1为真核特有", "both"),
    ("Rap1 signaling", "signaling", "Rap1通路为真核特有", "both"),
    ("Ras signaling", "signaling", "Ras通路为真核特有", "both"),
    ("Calcium signaling", "signaling", "钙信号为真核特有", "both"),
    ("cAMP signaling", "signaling", "真核cAMP信号（原核用cAMP-CRP不同）", "both"),
    ("neuroactive ligand-receptor", "signaling", "神经配体-受体互作为动物特有", "both"),
    ("cytokine-cytokine receptor", "signaling", "细胞因子-受体为动物特有", "both"),

    # === 6. 免疫系统 (Immune system) ===
    ("complement and coagulation", "immune", "补体与凝血级联为动物特有", "both"),
    ("Toll-like receptor", "immune", "TLR为动物先天免疫特有", "both"),
    ("NOD-like receptor", "immune", "NLR为动物先天免疫特有", "both"),
    ("RIG-I-like receptor", "immune", "RLR为动物先天免疫特有", "both"),
    ("T cell receptor", "immune", "TCR为脊椎动物适应性免疫特有", "both"),
    ("B cell receptor", "immune", "BCR为脊椎动物适应性免疫特有", "both"),
    ("antigen processing", "immune", "抗原加工呈递为动物特有", "both"),
    ("NK cell", "immune", "NK细胞为脊椎动物特有", "both"),
    ("Th1", "immune", "Th1极化为脊椎动物特有", "both"),
    ("Th2", "immune", "Th2极化为脊椎动物特有", "both"),
    ("Th17", "immune", "Th17极化为脊椎动物特有", "both"),

    # === 7. 神经系统 (Nervous system) ===
    ("synaptic vesicle", "nervous", "突触囊泡为神经元特有", "both"),
    ("glutamatergic synapse", "nervous", "谷氨酸能突触为动物特有", "both"),
    ("dopaminergic synapse", "nervous", "多巴胺能突触为动物特有", "both"),
    ("GABAergic synapse", "nervous", "GABA能突触为动物特有", "both"),
    ("cholinergic synapse", "nervous", "胆碱能突触为动物特有", "both"),
    ("serotonergic synapse", "nervous", "5-HT能突触为动物特有", "both"),
    ("long-term potentiation", "nervous", "LTP为神经元特有", "both"),
    ("long-term depression", "nervous", "LTD为神经元特有", "both"),

    # === 8. 内分泌 (Endocrine system) ===
    ("thyroid hormone", "endocrine", "甲状腺激素信号为脊椎动物特有", "both"),
    ("oxytocin", "endocrine", "催产素信号为脊椎动物特有", "both"),
    ("prolactin", "endocrine", "催乳素信号为脊椎动物特有", "both"),
    ("GnRH", "endocrine", "GnRH信号为脊椎动物特有", "both"),
    ("vasopressin", "endocrine", "加压素信号为脊椎动物特有", "both"),
    ("melanogenesis", "endocrine", "黑素生成为动物特有", "both"),
    ("adipocytokine", "endocrine", "脂肪因子信号为动物特有", "both"),
    ("progesterone", "endocrine", "孕酮信号为动物特有", "both"),
    ("estrogen", "endocrine", "雌激素信号为动物特有", "both"),
    ("androgen", "endocrine", "雄激素信号为动物特有", "both"),

    # === 9. 感觉系统 (Sensory systems) ===
    ("olfactory transduction", "sensory", "嗅觉转导为动物特有", "both"),
    ("taste transduction", "sensory", "味觉转导为动物特有", "both"),
    ("phototransduction", "sensory", "光转导（视杆/视锥）为动物特有", "both"),

    # === 10. 蛋白降解与自噬/死亡 (Proteolysis, autophagy, cell death) ===
    ("ubiquitin", "degradation_autophagy", "泛素化降解系统为真核特有", "both"),
    ("autophagy", "degradation_autophagy", "自噬为真核特有", "both"),
    ("apoptosis", "degradation_autophagy", "细胞凋亡为真核特有", "both"),
    ("necroptosis", "degradation_autophagy", "程序性坏死为动物特有", "both"),

    # === 11. 脂质代谢（真核特有分支） (Lipid metabolism) ===
    ("cholesterol metabolism", "lipid", "胆固醇代谢为真核特有", "both"),
    ("sphingolipid signaling", "lipid", "鞘脂信号为真核特有", "both"),

    # === 12. 疾病-神经退行 (Neurodegenerative disease) ===
    ("Alzheimer disease", "disease_neuro", "阿尔茨海默病为人类神经退行病", "both"),
    ("Parkinson disease", "disease_neuro", "帕金森病为人类神经退行病", "both"),
    ("Huntington disease", "disease_neuro", "亨廷顿病为人类神经退行病", "both"),
    ("amyotrophic lateral sclerosis", "disease_neuro", "ALS为人类神经退行病", "both"),
    ("pathways of neurodegeneration", "disease_neuro", "神经退行性疾病集合通路", "both"),

    # === 13. 疾病-感染 (Infectious disease, human-specific) ===
    ("human ", "disease_infection", "人类宿主特异性通路（含空格避免误伤）", "both"),
    ("coronavirus disease", "disease_infection", "COVID-19等冠状病毒病", "both"),
    ("influenza A", "disease_infection", "甲型流感", "both"),
    ("hepatitis ", "disease_infection", "肝炎（含空格避免误伤）", "module"),
    ("hepatitis C", "disease_infection", "丙型肝炎", "pathway"),
    ("measles", "disease_infection", "麻疹", "both"),
    ("malaria", "disease_infection", "疟疾", "both"),

    # === 14. 疾病-癌症 (Cancer) ===
    ("cancer:", "disease_cancer", "癌症泛通路（含冒号避免误伤）", "pathway"),
    ("colorectal cancer", "disease_cancer", "结直肠癌", "pathway"),
    ("breast cancer", "disease_cancer", "乳腺癌", "pathway"),
    ("prostate cancer", "disease_cancer", "前列腺癌", "pathway"),
    ("melanoma", "disease_cancer", "黑色素瘤", "pathway"),
    ("acute myeloid leukemia", "disease_cancer", "急性髓系白血病", "pathway"),
    ("chronic myeloid leukemia", "disease_cancer", "慢性髓系白血病", "pathway"),
    ("renal cell carcinoma", "disease_cancer", "肾细胞癌", "pathway"),
    ("pancreatic cancer", "disease_cancer", "胰腺癌", "pathway"),
    ("gastric cancer", "disease_cancer", "胃癌", "pathway"),
    ("hepatocellular carcinoma", "disease_cancer", "肝细胞癌", "pathway"),
    ("non-small cell lung cancer", "disease_cancer", "非小细胞肺癌", "pathway"),
    ("small cell lung cancer", "disease_cancer", "小细胞肺癌", "pathway"),
    ("thyroid cancer", "disease_cancer", "甲状腺癌", "pathway"),
    ("endometrial cancer", "disease_cancer", "子宫内膜癌", "pathway"),
    ("basal cell carcinoma", "disease_cancer", "基底细胞癌", "pathway"),
    ("viral carcinogenesis", "disease_cancer", "病毒致癌", "pathway"),
    ("chemical carcinogenesis", "disease_cancer", "化学致癌", "pathway"),

    # === 15. 疾病-代谢/免疫/心血管 (Metabolic / immune / CV disease) ===
    ("diabetes", "disease_metabolic", "糖尿病为人类代谢病", "both"),
    ("obesity", "disease_metabolic", "肥胖为人类代谢病", "both"),
    ("inflammatory bowel disease", "disease_metabolic", "IBD为人类免疫病", "both"),
    ("atherosclerosis", "disease_metabolic", "动脉粥样硬化为人类心血管病", "both"),
    ("asthma", "disease_metabolic", "哮喘为人类呼吸病", "both"),
]


# ---------------------------------------------------------------------------
# Pathway ID / prefix blacklists (not keyword-based)
# ---------------------------------------------------------------------------
# ko01100/01110/01120 are global metabolic maps (not real pathways).
PATHWAY_EXCLUDE_IDS: frozenset[str] = frozenset({
    "ko01100",
    "ko01110",
    "ko01120",
})

# ko041xx-049xx, ko05xxx blocks are mostly eukaryote cellular processes,
# organismal systems, and human diseases in KEGG pathway classification.
PATHWAY_EXCLUDE_PREFIXES: tuple[str, ...] = (
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


# ---------------------------------------------------------------------------
# Derived keyword tuples (backward-compatible with build_prok_allowlists.py)
# ---------------------------------------------------------------------------

def _keywords_for(target: str) -> tuple[str, ...]:
    """Return keywords that apply to ``target`` ('module' | 'pathway').

    'both' entries apply to both; entries with matching applies_to apply;
    entries with the other single target do NOT apply.
    """
    out: list[str] = []
    for kw, _cat, _reason, applies in _BLACKLIST_ENTRIES:
        if applies == "both" or applies == target:
            out.append(kw)
    # Deduplicate while preserving order (some keywords may repeat across
    # categories in future).
    seen: set[str] = set()
    deduped: list[str] = []
    for kw in out:
        if kw not in seen:
            seen.add(kw)
            deduped.append(kw)
    return tuple(deduped)


MODULE_EXCLUDE_KEYWORDS: tuple[str, ...] = _keywords_for("module")
PATHWAY_EXCLUDE_KEYWORDS: tuple[str, ...] = _keywords_for("pathway")


# ---------------------------------------------------------------------------
# Public accessors
# ---------------------------------------------------------------------------

def all_blacklist_entries() -> list[BlacklistEntry]:
    """Return all blacklist entries as dataclasses."""
    return [BlacklistEntry(kw, cat, reason, appl) for kw, cat, reason, appl in _BLACKLIST_ENTRIES]


def blacklist_entries_for(target: str) -> list[BlacklistEntry]:
    """Return entries applicable to ``target`` ('module' | 'pathway')."""
    return [
        e for e in all_blacklist_entries()
        if e.applies_to == "both" or e.applies_to == target
    ]


def blacklist_summary() -> dict:
    """Return a summary dict of the blacklist (counts per category/target)."""
    from collections import Counter
    entries = all_blacklist_entries()
    by_cat = Counter(e.category for e in entries)
    by_target = Counter(e.applies_to for e in entries)
    return {
        "blacklist_version": BLACKLIST_VERSION,
        "total_entries": len(entries),
        "module_keyword_count": len(MODULE_EXCLUDE_KEYWORDS),
        "pathway_keyword_count": len(PATHWAY_EXCLUDE_KEYWORDS),
        "pathway_excluded_ids": sorted(PATHWAY_EXCLUDE_IDS),
        "pathway_excluded_prefixes": list(PATHWAY_EXCLUDE_PREFIXES),
        "entries_by_category": dict(by_cat),
        "entries_by_applies_to": dict(by_target),
    }


# ---------------------------------------------------------------------------
# Blacklist file writers (mirror allowlist output style)
# ---------------------------------------------------------------------------

def write_blacklist_files(out_dir: Path, version: str = BLACKLIST_VERSION) -> dict:
    """Write blacklist TSV + TXT files mirroring allowlist output style.

    Outputs (in ``out_dir``):
      - kegg_module_blacklist_v{version}.tsv  (keyword, category, reason)
      - kegg_module_blacklist_v{version}.txt  (one keyword per line)
      - kegg_pathway_blacklist_v{version}.tsv (rule_type, value, category, reason)
      - kegg_pathway_blacklist_v{version}.txt (one rule per line, prefixed)
      - kegg_blacklist_meta_v{version}.json   (summary + file list)
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- MODULE blacklist (keyword-only) ---
    mod_entries = blacklist_entries_for("module")
    mod_tsv = out_dir / f"kegg_module_blacklist_v{version}.tsv"
    with mod_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["keyword", "category", "reason", "applies_to"])
        for e in mod_entries:
            w.writerow([e.keyword, e.category, e.reason, e.applies_to])

    mod_txt = out_dir / f"kegg_module_blacklist_v{version}.txt"
    with mod_txt.open("w", encoding="utf-8") as f:
        f.write(f"# kegg module blacklist v{version}\n")
        f.write("# One keyword per line; case-insensitive substring match on module name\n")
        f.write(f"# total: {len(mod_entries)} keywords\n")
        for e in mod_entries:
            f.write(f"{e.keyword}\n")

    # --- PATHWAY blacklist (keyword + id + prefix) ---
    path_entries = blacklist_entries_for("pathway")
    path_tsv = out_dir / f"kegg_pathway_blacklist_v{version}.tsv"
    with path_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["rule_type", "value", "category", "reason"])
        # keyword rules first
        for e in path_entries:
            w.writerow(["keyword", e.keyword, e.category, e.reason])
        # id rules
        for ko_id in sorted(PATHWAY_EXCLUDE_IDS):
            w.writerow(["id", ko_id, "global_map", "全局代谢图，非独立通路"])
        # prefix rules
        for prefix in PATHWAY_EXCLUDE_PREFIXES:
            w.writerow(["prefix", prefix, "pathway_block",
                        "KEGG通路分类块（细胞过程/生物体系统/人类疾病）真核为主"])

    path_txt = out_dir / f"kegg_pathway_blacklist_v{version}.txt"
    with path_txt.open("w", encoding="utf-8") as f:
        f.write(f"# kegg pathway blacklist v{version}\n")
        f.write("# rule_type:value per line; keyword=substring, id=exact, prefix=prefix match\n")
        f.write(f"# keyword rules: {len(path_entries)}\n")
        f.write(f"# id rules: {len(PATHWAY_EXCLUDE_IDS)}\n")
        f.write(f"# prefix rules: {len(PATHWAY_EXCLUDE_PREFIXES)}\n")
        for e in path_entries:
            f.write(f"keyword:{e.keyword}\n")
        for ko_id in sorted(PATHWAY_EXCLUDE_IDS):
            f.write(f"id:{ko_id}\n")
        for prefix in PATHWAY_EXCLUDE_PREFIXES:
            f.write(f"prefix:{prefix}\n")

    # --- Meta JSON ---
    meta = blacklist_summary()
    meta["files"] = {
        "module_tsv": mod_tsv.name,
        "module_txt": mod_txt.name,
        "pathway_tsv": path_tsv.name,
        "pathway_txt": path_txt.name,
    }
    meta_path = out_dir / f"kegg_blacklist_meta_v{version}.json"
    meta_path.write_text(
        json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return meta


# ---------------------------------------------------------------------------
# Allowlist source loaders (uniform source access for allowlists too)
# ---------------------------------------------------------------------------

def _allowlist_path(out_dir: Path, domain: str, kind: str, version: str, ext: str) -> Path:
    """Resolve allowlist file path.

    domain: 'bacteria' | 'archaea' | 'prok'
    kind:   'module' | 'pathway'
    ext:    'tsv' | 'txt'
    """
    if domain in ("bacteria", "archaea"):
        if kind == "module":
            return out_dir / f"kegg_module_allowlist_{domain}_v{version}.tsv"
        # pathway
        return out_dir / f"prokaryote_ko_pathway_allowlist_{domain}_v{version}.{ext}"
    # legacy prok union
    if kind == "module":
        return out_dir / f"kegg_module_allowlist_prok_v{version}.tsv"
    return out_dir / f"prokaryote_ko_pathway_allowlist_v{version}.{ext}"


def load_allowlist(
    out_dir: Path,
    domain: str,
    kind: str,
    version: str,
    id_column: str | None = None,
) -> set[str]:
    """Load allowlist IDs from file. Returns a set of IDs.

    For txt files, returns non-comment lines. For tsv files, returns values
    in ``id_column`` (defaults: module_id / pathway_id).
    """
    if id_column is None:
        id_column = "module_id" if kind == "module" else "pathway_id"

    # Prefer tsv for module, txt for pathway (matches existing outputs)
    if kind == "module":
        path = _allowlist_path(out_dir, domain, kind, version, "tsv")
    else:
        path = _allowlist_path(out_dir, domain, kind, version, "txt")
        if not path.exists():
            path = _allowlist_path(out_dir, domain, kind, version, "tsv")

    if not path.exists():
        return set()

    if path.suffix == ".txt":
        return {
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        }
    ids: set[str] = set()
    with path.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            val = row.get(id_column, "").strip()
            if val:
                ids.add(val)
    return ids


def load_blacklist(out_dir: Path, target: str, version: str = BLACKLIST_VERSION) -> set[str]:
    """Load blacklist keywords/IDs for ``target`` ('module' | 'pathway').

    For 'pathway' also includes id: and prefix: rules encoded as
    'id:ko01100' / 'prefix:ko041' strings.
    """
    if target == "module":
        path = out_dir / f"kegg_module_blacklist_v{version}.txt"
    else:
        path = out_dir / f"kegg_pathway_blacklist_v{version}.txt"
    if not path.exists():
        return set()
    rules: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        rules.add(line)
    return rules


__all__ = [
    "BLACKLIST_VERSION",
    "BlacklistEntry",
    "MODULE_EXCLUDE_KEYWORDS",
    "PATHWAY_EXCLUDE_KEYWORDS",
    "PATHWAY_EXCLUDE_IDS",
    "PATHWAY_EXCLUDE_PREFIXES",
    "all_blacklist_entries",
    "blacklist_entries_for",
    "blacklist_summary",
    "write_blacklist_files",
    "load_allowlist",
    "load_blacklist",
]
