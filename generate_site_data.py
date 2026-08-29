#!/usr/bin/env python3
"""Generate flat JSON data files for the SensASP GitHub Pages site."""

import json
import os

SRC = "senescence_sasp_database.json"
OUT = "docs/data"
os.makedirs(OUT, exist_ok=True)

with open(SRC) as f:
    db = json.load(f)

metadata = db["metadata"]
genes_raw = db["genes"]

# ── 1. genes.json — flat records for Tabulator ──────────────────────────────
flat = []
for g in genes_raw:
    ids  = g["identifiers"]
    prov = g["provenance"]
    con  = g["conservation"]
    exp  = g["expression"]
    itx  = g["interactions"]

    top_tissues = exp.get("gtex_top_tissues", [])
    top_tissue_str = ", ".join(
        f"{t['tissue'].replace('_', ' ')} ({t['median_tpm']:.1f})"
        for t in top_tissues[:3]
    )

    partners = itx.get("string_top_partners", [])
    partners_str = ", ".join(p["partner"] for p in partners[:5])

    sources_list = prov.get("sources", [])
    flat.append({
        "symbol":              ids["hgnc_symbol"],
        "gene_name":           ids.get("gene_name", ""),
        "ensembl":             ids.get("ensembl_gene", ""),
        "entrez":              ids.get("entrez", ""),
        "uniprot":             ids.get("uniprot", ""),
        "sources":             ", ".join(sources_list),
        "n_sources":           prov.get("n_sources", 0),
        "in_cellage":          "CellAge"  in sources_list,
        "in_genage":           "GenAge"   in sources_list,
        "in_senmayo":          "SenMayo"  in sources_list,
        "in_reactome":         "Reactome" in sources_list,
        "conserved_mouse":     con.get("conserved_in_mouse", False),
        "mouse_ortholog":      con.get("mouse_ortholog_symbol", ""),
        "mouse_identity":      round(con["mouse_perc_identity"], 1) if con.get("mouse_perc_identity") else None,
        "conserved_zebrafish": con.get("conserved_in_zebrafish", False),
        "gtex_top_tissue":     top_tissues[0]["tissue"].replace("_", " ") if top_tissues else "",
        "gtex_max_tpm":        exp.get("gtex_max_tpm"),
        "gtex_median_tpm":     exp.get("gtex_median_across_tissues"),
        "gtex_n_tissues":      exp.get("gtex_n_tissues_detected_tpm1"),
        "gtex_top3_tissues":   top_tissue_str,
        "hpa_tissue_spec":     exp.get("hpa_rna_tissue_specificity", ""),
        "hpa_cell_spec":       exp.get("hpa_rna_single_cell_specificity", ""),
        "string_degree":       itx.get("string_degree_within_set", 0),
        "string_top_partners": partners_str,
    })

with open(f"{OUT}/genes.json", "w") as f:
    json.dump(flat, f, separators=(",", ":"))
print(f"genes.json: {len(flat)} records")

# ── 2. metrics.json — summary stats ─────────────────────────────────────────
with open("senescence_db_pipeline/results/build_metrics.json") as f:
    metrics = json.load(f)

out_metrics = {
    "n_genes": metadata["n_genes"],
    "build_date": metadata["build_date"],
    "sources": metrics["sources"],
    "n_sources_distribution": metrics["n_sources_distribution"],
    "multi_source_genes": metrics["multi_source_genes"],
    "coverage": metrics["coverage"],
}
with open(f"{OUT}/metrics.json", "w") as f:
    json.dump(out_metrics, f, indent=2)
print("metrics.json written")

# ── 3. sources.json — seed source metadata ───────────────────────────────────
sources = []
for name, info in metadata["seed_sources"].items():
    sources.append({
        "name":        name,
        "n_genes":     info["n_input"],
        "url":         info["url"],
        "description": info["description"],
    })
with open(f"{OUT}/sources.json", "w") as f:
    json.dump(sources, f, indent=2)
print("sources.json written")

# ── 4. genes_full.json — full nested records for gene detail pages ───────
with open(f"{OUT}/genes_full.json", "w") as f:
    json.dump(genes_raw, f, separators=(",", ":"))
print(f"genes_full.json: {len(genes_raw)} records")

print("\nDone — all data files written to docs/data/")
