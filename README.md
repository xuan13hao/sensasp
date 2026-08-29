# Unified Cellular Senescence / SASP Gene Database — Pipeline

A fully scripted, one-command-reproducible pipeline that unifies four senescence /
SASP gene resources onto a single canonical identifier and enriches every gene
with cross-species conservation, tissue/cell-type expression, and protein–protein
interaction annotation.

## Quick start

```bash
./main.sh            # install deps, run the data pipeline, regenerate all figures
./main.sh --no-figs  # data pipeline only (skip the notebook)
```

Everything is fetched from public APIs at run time — no manual downloads, no
credentials. Total runtime is a few minutes (the Human Protein Atlas step makes
one request per gene and dominates the wall time).

## What it does

| Step | Module | Output |
|---|---|---|
| 1 | `step1_acquire.py` | Download CellAge, GenAge, SenMayo, Reactome → `data/raw/seed_lists_raw.json` |
| 2 | `step2_harmonize.py` | Map symbols→Ensembl/UniProt/Entrez via MyGene.info → `data/processed/id_map.{json,csv}` |
| 3 | `step3_unify.py` | Collapse to unique genes on Ensembl ID, record provenance (RQ1) → `data/processed/unified_core.json`, `results/source_overlap.csv` |
| 4 | `step4_conservation.py` | Mouse + zebrafish orthologs via Ensembl BioMart → `data/processed/conservation.json` |
| 5 | `step5_expression.py` | GTEx v8 (54 tissues) + Human Protein Atlas → `data/processed/expression.json` |
| 6 | `step6_interactions.py` | STRING v12 high-confidence network (score≥700) → `data/processed/interaction.json`, `string_network.json` |
| 7 | `step7_assemble.py` | Merge all layers → `results/senescence_sasp_database.{json,csv}` |
| 8 | `step8_metrics.py` | Outcome metrics (RQ1/RQ2) → `results/build_metrics.json`, `summary_report.md` |

Figures are generated separately from the pipeline outputs by
`notebooks/figures.ipynb`.

## Running individual steps

```bash
PYTHONPATH=. python run_pipeline.py          # all steps
PYTHONPATH=. python run_pipeline.py 4 5 6     # only these (reuse earlier outputs)
```

## Layout

```
senescence_db_pipeline/
├── main.sh                 # one-click reproduction
├── run_pipeline.py         # step runner
├── requirements.txt
├── pipeline/               # step1..step8 modules + config.py
├── notebooks/figures.ipynb # all figure-generation code
├── data/raw, data/processed
├── results/                # database + metrics + report
├── figures/                # generated PNGs
└── paper/paper.md          # manuscript (Intro/Methods/Data/Results/Discussion)
```

## Data sources

CellAge & GenAge (HAGR), SenMayo (MSigDB `SAUL_SEN_MAYO` / M45803, Saul et al.
2022), Reactome Cellular Senescence (R-HSA-2559583); annotation via MyGene.info,
Ensembl BioMart, GTEx v8, Human Protein Atlas, and STRING v12.

## Web database (GitHub Pages)

The interactive database is served from `docs/` via GitHub Pages:

```
https://xuan13hao.github.io/sensasp
```

To enable: **Settings → Pages → Source: Deploy from branch → Branch: `main` → Folder: `/docs`**

To regenerate site data after re-running the pipeline:

```bash
python3 generate_site_data.py
```

This rebuilds `docs/data/genes.json`, `genes_full.json`, `metrics.json`, and `sources.json`.

## Key results

1,460 summed source entries → **1,250 unique genes** (210 redundant collapsed,
14.4%); 173 genes confirmed by ≥2 sources, 2 (IL6, JUN) by all four. Net-new
annotation reaches 95.8% (conservation), 97.8% (expression), 92.9%
(interactions), with 89.3% complete across all three layers.
