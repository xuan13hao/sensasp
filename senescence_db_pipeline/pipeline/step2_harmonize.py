"""Step 2 — Harmonize gene identifiers.

Maps every seed symbol (with alias resolution) to a canonical identifier triple
(HGNC symbol, Ensembl gene, UniProt) via MyGene.info, and writes id_map.json/csv.
"""
import json
import csv
from collections import defaultdict

import requests

from . import config as C


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def query_mygene(symbols):
    """POST symbols to MyGene.info in chunks; return concatenated record list."""
    fields = "symbol,name,taxid,entrezgene,ensembl.gene,uniprot,HGNC"
    records = []
    for ch in _chunks(symbols, C.MYGENE_CHUNK):
        resp = requests.post(
            C.URL_MYGENE,
            json={"q": ch, "scopes": "symbol,alias", "fields": fields, "species": "human"},
            timeout=C.REQUEST_TIMEOUT,
        )
        resp.raise_for_status()
        records.extend(resp.json())
    return records


def _get_ensembl(r):
    e = r.get("ensembl")
    if isinstance(e, list):
        return e[0].get("gene") if e else None
    if isinstance(e, dict):
        return e.get("gene")
    return None


def _get_uniprot(r):
    u = r.get("uniprot")
    if isinstance(u, dict):
        sp = u.get("Swiss-Prot")
        return sp[0] if isinstance(sp, list) and sp else (sp if isinstance(sp, str) else None)
    return None


def _pick(records, qsym):
    """Choose the best record for an input symbol: human, exact-symbol, has Ensembl, top score."""
    human = [r for r in records if r.get("taxid") == 9606]
    pool = human or records
    exact = [r for r in pool if r.get("symbol", "").upper() == qsym.upper()]
    cand = exact or pool
    return sorted(cand, key=lambda r: (0 if r.get("ensembl") else 1, -(r.get("_score") or 0)))[0]


def run():
    seed = json.loads((C.DATA_RAW / "seed_lists_raw.json").read_text())
    all_syms = sorted(set().union(*[set(v["genes"]) for v in seed.values()]))

    records = query_mygene(all_syms)
    by_q = defaultdict(list)
    for r in records:
        if "query" in r:
            by_q[r["query"]].append(r)

    id_map = {}
    for qsym in all_syms:
        recs = by_q.get(qsym, [])
        if not recs or all(r.get("notfound") for r in recs):
            id_map[qsym] = {"input_symbol": qsym, "hgnc_symbol": None, "gene_name": None,
                            "entrez": None, "ensembl_gene": None, "uniprot": None,
                            "hgnc_id": None, "symbol_changed": False, "mapped": False}
            continue
        r = _pick(recs, qsym)
        hgnc = r.get("HGNC")
        id_map[qsym] = {
            "input_symbol": qsym,
            "hgnc_symbol": r.get("symbol"),
            "gene_name": r.get("name"),
            "entrez": r.get("entrezgene"),
            "ensembl_gene": _get_ensembl(r),
            "uniprot": _get_uniprot(r),
            "hgnc_id": f"HGNC:{hgnc}" if hgnc else None,
            "symbol_changed": r.get("symbol", "").upper() != qsym.upper(),
            "mapped": True,
        }

    (C.DATA_PROC / "id_map.json").write_text(json.dumps(id_map, indent=2))
    with open(C.DATA_PROC / "id_map.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(next(iter(id_map.values())).keys()))
        w.writeheader()
        for v in sorted(id_map.values(), key=lambda x: x["input_symbol"]):
            w.writerow(v)

    n_ens = sum(1 for v in id_map.values() if v["ensembl_gene"])
    n_uni = sum(1 for v in id_map.values() if v["uniprot"])
    n_chg = sum(1 for v in id_map.values() if v["symbol_changed"])
    print(f"[step2] {len(id_map)} symbols | Ensembl={n_ens} UniProt={n_uni} symbol_updated={n_chg}")
    return id_map


if __name__ == "__main__":
    run()
