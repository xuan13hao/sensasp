"""Step 5 — Add tissue/cell-type expression (GTEx v8 + Human Protein Atlas).

Resolves Ensembl gene ids to versioned GENCODE ids, pulls GTEx median TPM across
all 54 tissues, and adds HPA RNA tissue/single-cell specificity + subcellular
location. Writes expression.json.
"""
import json
import time
from collections import defaultdict

import requests

from . import config as C


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def resolve_gencode(ens_ids):
    """Map ENSG -> versioned GENCODE id (gencodeId) via the GTEx reference/gene endpoint."""
    mapping = {}
    for ch in _chunks(ens_ids, 200):
        r = requests.get(f"{C.URL_GTEX}/reference/gene",
                         params={"geneId": ch, "datasetId": C.GTEX_DATASET,
                                 "page": 0, "itemsPerPage": 100000},
                         timeout=C.REQUEST_TIMEOUT)
        r.raise_for_status()
        for rec in r.json().get("data", []):
            g = rec.get("gencodeId")
            base = (rec.get("geneId") or (g.split(".")[0] if g else None))
            if base and g:
                mapping[base] = g
        time.sleep(0.2)
    return mapping


def fetch_gtex_medians(gencode_ids):
    """Median TPM per tissue for each GENCODE id."""
    out = defaultdict(dict)
    for ch in _chunks(gencode_ids, 200):
        r = requests.get(f"{C.URL_GTEX}/expression/medianGeneExpression",
                         params={"gencodeId": ch, "datasetId": C.GTEX_DATASET,
                                 "page": 0, "itemsPerPage": 100000},
                         timeout=C.REQUEST_TIMEOUT)
        r.raise_for_status()
        for rec in r.json().get("data", []):
            out[rec["gencodeId"]][rec["tissueSiteDetailId"]] = rec["median"]
        time.sleep(0.2)
    return out


def fetch_hpa(ens_ids):
    """HPA RNA specificity + subcellular location, one gene per call (search_download)."""
    cols = "g,eg,rnats,rnatd,rnascs,scl"
    out = {}
    for i, e in enumerate(ens_ids):
        try:
            r = requests.get(C.URL_HPA,
                             params={"search": e, "format": "json",
                                     "columns": cols, "compress": "no"},
                             timeout=C.REQUEST_TIMEOUT)
            r.raise_for_status()
            rows = r.json()
            rec = next((x for x in rows if x.get("Ensembl") == e), rows[0] if rows else None)
            if rec:
                out[e] = {
                    "rna_tissue_specificity": rec.get("RNA tissue specificity"),
                    "rna_tissue_distribution": rec.get("RNA tissue distribution"),
                    "rna_singlecell_specificity": rec.get("RNA single cell type specificity"),
                    "subcellular_location": rec.get("Subcellular location"),
                }
        except Exception:
            pass
        time.sleep(C.REQUEST_PAUSE)
        if (i + 1) % 200 == 0:
            print(f"[step5]   HPA {i+1}/{len(ens_ids)}")
    return out


def run():
    unified = json.loads((C.DATA_PROC / "unified_core.json").read_text())
    ens_ids = sorted({g["ensembl_gene"] for g in unified.values() if g.get("ensembl_gene")})

    gencode = resolve_gencode(ens_ids)
    medians = fetch_gtex_medians(sorted(set(gencode.values())))
    hpa = fetch_hpa(ens_ids)

    expr = {}
    for gid, g in unified.items():
        e = g.get("ensembl_gene")
        gc = gencode.get(e)
        expr[gid] = {
            "gencode_id": gc,
            "gtex_median_tpm": medians.get(gc) if gc else None,
            "hpa": hpa.get(e),
        }

    (C.DATA_PROC / "expression.json").write_text(json.dumps(expr, indent=2))
    n_gtex = sum(1 for v in expr.values() if v["gtex_median_tpm"])
    n_hpa = sum(1 for v in expr.values() if v["hpa"])
    n = len(expr)
    print(f"[step5] GTEx={n_gtex} ({100*n_gtex/n:.1f}%) HPA={n_hpa} ({100*n_hpa/n:.1f}%)")
    return expr


if __name__ == "__main__":
    run()
