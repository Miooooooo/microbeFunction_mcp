# mgnify_mcp

MGnify 物种索引检索与注释按需下载的 MCP 服务（stdio）。

## 数据路径

默认读取仓库根目录下：

- `data/mgnify/release_manifest.json`
- `data/mgnify/{biome}_{release}_species_index.jsonl`
- 注释缓存：`downloads/{species_rep}/genome/`

可通过环境变量 `MGNIFY_DATA_DIR` 覆盖索引目录（仍使用仓库根目录解析 `cache_dir`）。

构建索引：

```bash
python build_mgnify_index.py --biome human-gut --release v2.0
```

## 启动 MCP（stdio）

在仓库根目录执行：

```bash
python -m mgnify_mcp.server
```

## Cursor 配置

本仓库已包含项目级配置：`.cursor/mcp.json`（打开 Pathway 工作区后自动加载）。

1. **Cursor Settings** → **MCP**，确认出现 `mgnify_mcp` 且状态为已连接（绿点）。
2. 若未出现：命令面板执行 **Developer: Reload Window**，或关闭后重新打开本仓库。
3. 本机 Python 路径与 `kegg_mcp` 一致时，无需改 `command`；若启动失败，把 `.cursor/mcp.json` 里的 `command` 改成你本机 `python.exe` 的绝对路径（与 `~/.cursor/mcp.json` 中 `kegg_mcp` 相同即可）。

可选：合并到用户全局 `~/.cursor/mcp.json` 的 `mcpServers` 中（在其他仓库也想用时）：

```json
"mgnify_mcp": {
  "command": "D:/Program Files/python.exe",
  "args": ["-m", "mgnify_mcp.server"],
  "cwd": "E:/LLM_learning/Project/Pathway",
  "env": {
    "PYTHONPATH": "E:/LLM_learning/Project/Pathway",
    "MGNIFY_DATA_DIR": "E:/LLM_learning/Project/Pathway/data/mgnify"
  }
}
```

**与 kegg_mcp 联用**：全局或项目 MCP 中同时启用 `kegg_mcp` 与 `mgnify_mcp` 即可；典型流程为 `search_species` → `fetch_annotations` → `kegg_module_completeness`。

## Tools

- `search_species`：`biome`, `species_query`, `release`（默认 `latest`）, `limit`
- `fetch_annotations`：`species_rep`, `biome`, `release`, `roles`, `convert_gtf`, `preview_rows`

## CLI

```bash
python -m mgnify_mcp search --biome human-gut --query "Enterobacter kobei"
python -m mgnify_mcp fetch --species-rep MGYG000000238 --biome human-gut --roles eggnog_tsv
```

## Resource

- `resource://mgnify/releases`：release_manifest 摘要（不含 JSONL 全文）
