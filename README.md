# microbeFunction_mcp

微生物功能分析 MCP 工具包

支持[mgnify_genomes](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/)数据库所有菌种(包含人肠道，口腔，阴道；小鼠肠道；海洋淤泥；土壤；根际；牛羊瘤胃等)的基因组信息下载和查询。

支持KEGG所有数据库查询和功能注释，代谢模块完整度分析。并构建细菌和古菌KEGG pathway/module的白名单和黑名单(只在真核中出现)及关键词黑名单用于进一步过滤在原核生物中不可能出现的功能模块。

支持COG数据库中的功能查询。

## 概览

| 工具包 | 功能 | 端口 |
|--------|------|------|
| `tools/kegg` | KEGG REST API + 代谢模块完整度分析 + COG数据库（依赖 `kegg-pathways-completeness` CLI） | 8791 |
| `tools/mgnify` | MGnify 物种索引检索 + 注释文件按需下载 | 8792 |
| `deploy/` | 统一 FastAPI 部署层，同时挂载两个 MCP | 8788 |

## 部署

### 1. 安装依赖

```bash
pip install uv
uv venv .venv --python=3.13
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
uv sync
```

### 2. 配置

```bash
cp default.conf.toml local.conf.toml
# 按需修改 local.conf.toml 中的端口等配置
```

### 3. 启动统一服务

```bash
uv run -m deploy.web
```

启动后：
- KEGG MCP: `http://127.0.0.1:8788/kegg_mcp/mcp/`
- MGnify MCP: `http://127.0.0.1:8788/mgnify_mcp/mcp/`
- 服务列表: `http://127.0.0.1:8788/api/list_mcps`

### 4. 单独启动某个 MCP

```bash
# KEGG (端口 8791)
uv run -m tools.kegg.deploy

# MGnify (端口 8792)
uv run -m tools.mgnify.deploy
```

## 典型工作流

```
search_species (mgnify)  →  fetch_annotations (mgnify)  →  kegg_module_completeness (kegg)   →   functional profiling of bacterial genomes
     找到物种 MAG              下载 eggNOG/GFF 注释            分析代谢模块完整度                         菌种基因组的功能画像
```

## 使用示例：生成 *Clostridium baratii* 功能画像

### 步骤 1：搜索物种

```bash
uv run -m tools.mgnify search --biome human-gut --query "Clostridium baratii"
```

输出示例：
```json
{
  "items": [
    {
      "species_rep": "MGYG000000064",
      "species_name": "Clostridium baratii",
      "completeness": 99.19,
      "contamination": 1.61,
      "genome_count": 12
    }
  ]
}
```

### 步骤 2：下载注释

```bash
uv run -m tools.mgnify fetch --species-rep MGYG000000064 --biome human-gut --roles eggnog_tsv,gff
```

注释文件会保存在 `downloads/MGYG000000064/genome/` 目录下。

### 步骤 3：分析代谢模块完整度

```bash
uv run -m tools.kegg analyze \
  --annotation-file downloads/MGYG000000064/genome/MGYG000000064_eggNOG.tsv \
  --kegg-column KEGG_ko \
  --output downloads/MGYG000000064/genome/MGYG000000064_module_completeness.tsv
```

输出示例：
```
output=downloads/MGYG000000064/genome/MGYG000000064_module_completeness.tsv 
unique_ko_count=1384 modules_with_any_hit=180 modules_above_threshold=180
```

### 步骤 4：解读结果

| 功能分类 | 完整模块数/总模块数 |
|----------|---------------------|
| Carbohydrate metabolism | 8/30 |
| Amino acid metabolism | 12/39 |
| Energy metabolism | 4/23 |
| Nucleotide metabolism | 6/10 |
| Glycan metabolism | 6/15 |
| Cofactor and vitamin metabolism | 7/38 |
| Lipid metabolism | 3/9 |

**代表性完整模块：**
- M00001: 糖酵解 (Embden-Meyerhof pathway)
- M00632: 半乳糖降解 (Leloir pathway)
- M00579: 乙酸生成 (磷酸乙酰转移酶-乙酸激酶)
- M00651: 万古霉素耐药 (D-Ala-D-Lac type)
- M00122: 钴胺素生物合成 (维生素B12)
- M00924: 钴胺素生物合成 (厌氧途径)

**部分完整模块：**
- M00003: 糖酵解 (78.57%)
- M00010: 乙醇发酵 (66.67%)
- M00009: 乳酸发酵 (50.0%)

**基因组概览：**
- 总基因数：2,893
- 有 KO 注释的基因：1,647 (56%)
- COG 分类：碳水化合物代谢 (229)、氨基酸代谢 (191)、能量代谢 (171)
- CAZy 酶：糖基转移酶 (GT) 22 个、糖苷水解酶 (GH) 16 个

## CLI 用法

```bash
# KEGG 模块完整度分析
uv run -m tools.kegg analyze --annotation-file downloads/MGYG000000238/genome/MGYG000000238_gene_annotations.tsv

# KEGG 批处理
uv run -m tools.kegg batch --manifest manifest.tsv --jobs 4

# MGnify 物种搜索
uv run -m tools.mgnify search --biome human-gut --query "Enterobacter kobei"

# MGnify 注释下载
uv run -m tools.mgnify fetch --species-rep MGYG000000238 --biome human-gut --roles eggnog_tsv
```

## 数据路径

- MGnify 索引: `data/mgnify/`
- MGnify 注释缓存: `downloads/{species_rep}/genome/`
- KEGG allowlist 数据: `tools/kegg/data/`
- COG 目录: `tools/kegg/data/COG.csv`

可通过环境变量 `MGNIFY_DATA_DIR` 覆盖 MGnify 索引目录。


## Cursor / Claude Desktop 配置

```json
{
  "mcpServers": {
    "kegg_mcp": {
      "url": "http://127.0.0.1:8788/kegg_mcp/mcp/",
      "transport": "streamable_http"
    },
    "mgnify_mcp": {
      "url": "http://127.0.0.1:8788/mgnify_mcp/mcp/",
      "transport": "streamable_http"
    }
  }
}
```
