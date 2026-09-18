# SenSASP: A Unified, Multi-Layer Database of Senescence and SASP Genes

[![PyPI](https://img.shields.io/pypi/v/sensasp-mcp)](https://pypi.org/project/sensasp-mcp/)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.08.26.747427-b31b1b)](https://doi.org/10.64898/2026.08.26.747427)
[![GitHub Pages](https://img.shields.io/badge/Database-xuan13hao.github.io%2Fsensasp-1d6fa4)](https://xuan13hao.github.io/sensasp)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)

**1,250 unique human senescence / SASP genes** unified from CellAge, GenAge, SenMayo, and
Reactome, enriched with cross-species conservation, GTEx expression, and STRING interactions.

> Xuan, H., Huang, Y., & Bian, J. (2026). SenSASP: A Unified, Multi-Layer Database of
> Senescence and SASP Genes. *bioRxiv* 2026.08.26.747427.
> https://doi.org/10.64898/2026.08.26.747427

---

## Web database

Browse and search all genes interactively — no account or software required:

**https://xuan13hao.github.io/sensasp**

- Searchable/sortable/filterable gene table (Tabulator.js)
- Per-gene detail pages: `sensasp/gene.html?gene=TP53`
- Coverage and source statistics
- One-click citation export (APA + BibTeX)

---

## MCP server (Claude AI integration)

SenSASP ships as an [MCP](https://modelcontextprotocol.io) server so Claude and other
AI assistants can query the database in natural language.

### Install

```bash
pip install sensasp-mcp
```

### Connect to Claude Desktop

Add to your Claude Desktop config file:

- **Linux:** `~/.config/claude/claude_desktop_config.json`
- **Mac:** `~/Library/Application Support/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "sensasp": {
      "command": "python3",
      "args": ["-m", "sensasp_mcp.server"]
    }
  }
}
```

Restart Claude Desktop. Then ask in natural language:

> *"Which SenSASP genes are confirmed by all four sources?"*
> *"Show me the full annotation for CDKN2A."*
> *"What are the top hub genes by STRING interaction degree?"*

### Zero-install with uvx

```bash
uvx sensasp-mcp
```

### Available tools

| Tool | Description |
|------|-------------|
| `search_genes` | Free-text + filter search (source, conservation, expression, STRING degree) |
| `get_gene` | Full record for one gene: conservation, GTEx tissues, STRING partners |
| `get_database_stats` | Gene count, coverage percentages, build metadata |
| `list_sources` | The 4 seed databases with gene counts and URLs |
| `get_top_hub_genes` | Top-N genes by STRING intra-set interaction degree |
| `get_multi_source_genes` | All 173 genes confirmed by ≥2 sources |

### Remote / self-hosted

Run as an HTTP server that any MCP client can reach without a local install:

```bash
python3 -m sensasp_mcp.server --transport sse --port 8080
```

Users connect via:

```json
{ "mcpServers": { "sensasp": { "url": "https://your-server.com/sse" } } }
```

---

## Pipeline

A fully scripted pipeline that reproduces the database from scratch:

```bash
./main.sh            # install deps, run all steps, regenerate figures
./main.sh --no-figs  # data only
```

Everything is fetched from public APIs at run time — no manual downloads.

### Steps

| Step | Module | Output |
|------|--------|--------|
| 1 | `step1_acquire.py` | CellAge, GenAge, SenMayo, Reactome → `data/raw/` |
| 2 | `step2_harmonize.py` | Symbol → Ensembl/UniProt/Entrez via MyGene.info |
| 3 | `step3_unify.py` | Collapse to unique genes, record provenance |
| 4 | `step4_conservation.py` | Mouse + zebrafish orthologs (Ensembl BioMart) |
| 5 | `step5_expression.py` | GTEx v8 (54 tissues) + Human Protein Atlas |
| 6 | `step6_interactions.py` | STRING v12 high-confidence PPI (score ≥ 700) |
| 7 | `step7_assemble.py` | Merge all layers → `senescence_sasp_database.{json,csv}` |
| 8 | `step8_metrics.py` | Coverage metrics → `build_metrics.json` |

```bash
PYTHONPATH=. python3 run_pipeline.py          # all steps
PYTHONPATH=. python3 run_pipeline.py 4 5 6   # specific steps only
```

### Regenerate site data

After re-running the pipeline:

```bash
python3 generate_site_data.py
```

Rebuilds `docs/data/genes.json`, `genes_full.json`, `metrics.json`, `sources.json`,
and the bundled data in `sensasp_mcp/data/`.

---

## Repository layout

```
sensasp/
├── senescence_sasp_database.{json,csv}  # final database
├── generate_site_data.py                # regenerate docs/data/ + sensasp_mcp/data/
├── mcp_server.py                        # standalone MCP server (dev)
├── docs/                                # GitHub Pages site
│   ├── index.html                       # gene table
│   ├── gene.html                        # per-gene detail page
│   ├── cite.html                        # citation page
│   └── data/                            # JSON files served to the browser
├── sensasp_mcp/                         # PyPI package (pip install sensasp-mcp)
│   ├── server.py                        # MCP server with stdio + SSE transport
│   └── data/                            # bundled JSON for offline use
└── senescence_db_pipeline/              # reproducible pipeline
    ├── pipeline/                        # step1..step8 modules
    ├── notebooks/figures.ipynb
    └── results/                         # build metrics, summary report
```

---

## Key results

1,460 summed source entries → **1,250 unique genes** (210 redundant collapsed, 14.4%);
173 genes confirmed by ≥2 sources, 2 (IL6, JUN) by all four. Net-new annotation
reaches 95.8% (conservation), 97.9% (expression), 93.0% (interactions),
with 89.4% complete across all three layers.

---

## Citation

```bibtex
@article{Xuan2026.08.26.747427,
  author       = {Xuan, Hao and Huang, Yu and Bian, Jiang},
  title        = {SenSASP: A Unified, Multi-Layer Database of Senescence and SASP Genes},
  elocation-id = {2026.08.26.747427},
  year         = {2026},
  doi          = {10.64898/2026.08.26.747427},
  publisher    = {Cold Spring Harbor Laboratory},
  URL          = {https://www.biorxiv.org/content/early/2026/09/11/2026.08.26.747427},
  journal      = {bioRxiv}
}
```
