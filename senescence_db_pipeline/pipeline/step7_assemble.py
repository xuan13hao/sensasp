"""Step 7 — Assemble the unified JSON database.

Merges the unified core with the conservation, expression and interaction layers
into a single senescence_sasp_database.json (+ a flat CSV view).
"""
import json
import csv

from . import config as C


def run():
    unified = json.loads((C.DATA_PROC / "unified_core.json").read_text())
    cons = json.loads((C.DATA_PROC / "conservation.json").read_text())
    expr = json.loads((C.DATA_PROC / "expression.json").read_text())
    inter = json.loads((C.DATA_PROC / "interaction.json").read_text())

    db = {
        "metadata": {
            "title": "Unified Cellular Senescence / SASP Gene Database",
            "sources": C.SOURCES,
            "canonical_key": "Ensembl gene id",
            "annotation_layers": ["conservation", "expression", "interaction"],
            "conservation_species": ["mouse (Mus musculus)", "zebrafish (Danio rerio)"],
            "expression_sources": ["GTEx v8", "Human Protein Atlas"],
            "interaction_source": f"STRING v12 (score>={C.STRING_MIN_SCORE})",
        },
        "genes": {},
    }

    for gid, g in unified.items():
        db["genes"][gid] = {
            "canonical_id": g["canonical_id"],
            "symbol": g["symbol"],
            "gene_name": g.get("gene_name"),
            "identifiers": {
                "ensembl_gene": g.get("ensembl_gene"),
                "uniprot": g.get("uniprot"),
                "entrez": g.get("entrez"),
                "hgnc_id": g.get("hgnc_id"),
                "input_symbols": g.get("input_symbols"),
            },
            "sources": g["sources"],
            "n_sources": g["n_sources"],
            "conservation": cons.get(gid),
            "expression": expr.get(gid),
            "interaction": inter.get(gid),
        }

    out = C.RESULTS / "senescence_sasp_database.json"
    out.write_text(json.dumps(db, indent=2))

    # flat CSV view
    with open(C.RESULTS / "senescence_sasp_database.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["canonical_id", "symbol", "gene_name", "ensembl_gene", "uniprot",
                    "entrez", "hgnc_id", "sources", "n_sources",
                    "mouse_ortholog", "zebrafish_ortholog",
                    "gtex_covered", "hpa_specificity", "string_degree"])
        for gid, g in db["genes"].items():
            c = g["conservation"] or {}
            e = g["expression"] or {}
            hpa = (e.get("hpa") or {}) if e else {}
            it = g["interaction"] or {}
            w.writerow([
                g["canonical_id"], g["symbol"], g.get("gene_name"),
                g["identifiers"]["ensembl_gene"], g["identifiers"]["uniprot"],
                g["identifiers"]["entrez"], g["identifiers"]["hgnc_id"],
                "|".join(g["sources"]), g["n_sources"],
                bool(c.get("mouse")), bool(c.get("zebrafish")),
                bool(e.get("gtex_median_tpm")), hpa.get("rna_tissue_specificity"),
                (it.get("string_degree") if it else None),
            ])

    print(f"[step7] assembled {len(db['genes'])} genes -> {out}")
    return db


if __name__ == "__main__":
    run()
