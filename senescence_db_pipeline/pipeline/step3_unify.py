"""Step 3 — Unify and deduplicate (RQ1).

Collapses the four seed lists onto the canonical Ensembl gene ID, records the
source provenance of every gene, and writes the unified core plus overlap stats.
"""
import json
import csv
from collections import defaultdict

from . import config as C


def run():
    seed = json.loads((C.DATA_RAW / "seed_lists_raw.json").read_text())
    id_map = json.loads((C.DATA_PROC / "id_map.json").read_text())

    # canonical key: Ensembl gene id when available, else HGNC symbol, else input symbol
    def canon(sym):
        m = id_map.get(sym, {})
        return m.get("ensembl_gene") or m.get("hgnc_symbol") or sym

    genes = {}
    summed = 0
    for src in C.SOURCES:
        for sym in seed[src]["genes"]:
            summed += 1
            key = canon(sym)
            m = id_map.get(sym, {})
            if key not in genes:
                genes[key] = {
                    "canonical_id": key,
                    "symbol": m.get("hgnc_symbol") or sym,
                    "input_symbols": set(),
                    "ensembl_gene": m.get("ensembl_gene"),
                    "uniprot": m.get("uniprot"),
                    "entrez": m.get("entrez"),
                    "hgnc_id": m.get("hgnc_id"),
                    "gene_name": m.get("gene_name"),
                    "sources": set(),
                }
            genes[key]["sources"].add(src)
            genes[key]["input_symbols"].add(sym)

    for g in genes.values():
        g["sources"] = sorted(g["sources"])
        g["input_symbols"] = sorted(g["input_symbols"])
        g["n_sources"] = len(g["sources"])

    unified = {k: genes[k] for k in sorted(genes)}
    (C.DATA_PROC / "unified_core.json").write_text(json.dumps(unified, indent=2))

    # overlap table: source-membership matrix
    with open(C.RESULTS / "source_overlap.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["canonical_id", "symbol", "n_sources"] + C.SOURCES)
        for g in unified.values():
            w.writerow([g["canonical_id"], g["symbol"], g["n_sources"]]
                       + [int(s in g["sources"]) for s in C.SOURCES])

    n = len(unified)
    dist = defaultdict(int)
    for g in unified.values():
        dist[g["n_sources"]] += 1
    redundant = summed - n
    print(f"[step3] summed entries={summed} -> unique genes={n} "
          f"| collapsed {redundant} ({100*redundant/summed:.1f}%)")
    print(f"[step3] n_sources distribution: "
          + " ".join(f"{k}src={dist[k]}" for k in sorted(dist)))
    return unified


if __name__ == "__main__":
    run()
