"""Step 4 — Add cross-species conservation (mouse, zebrafish) via Ensembl BioMart.

Queries orthologs in chunks using the BioMart XML/martservice REST endpoint and
writes conservation.json keyed by canonical Ensembl gene id.
"""
import json
import io
import csv
import time

import requests

from . import config as C

_ORTHO = {
    "mouse": ("mmusculus", "mmusculus_homolog_ensembl_gene",
              "mmusculus_homolog_associated_gene_name",
              "mmusculus_homolog_orthology_type", "mmusculus_homolog_perc_id"),
    "zebrafish": ("drerio", "drerio_homolog_ensembl_gene",
                  "drerio_homolog_associated_gene_name",
                  "drerio_homolog_orthology_type", "drerio_homolog_perc_id"),
}


def _biomart_xml(ids, attrs):
    filt = ",".join(ids)
    attr_xml = "".join(f'<Attribute name="{a}"/>' for a in attrs)
    return (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="0" '
        'uniqueRows="1" count="" datasetConfigVersion="0.6">'
        '<Dataset name="hsapiens_gene_ensembl" interface="default">'
        f'<Filter name="ensembl_gene_id" value="{filt}"/>'
        '<Attribute name="ensembl_gene_id"/>'
        f'{attr_xml}</Dataset></Query>'
    )


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def fetch_orthologs(ens_ids, species):
    _, e_gene, e_name, e_type, e_perc = _ORTHO[species]
    attrs = [e_gene, e_name, e_type, e_perc]
    out = {}
    for ch in _chunks(ens_ids, C.BIOMART_CHUNK):
        xml = _biomart_xml(ch, attrs)
        r = requests.get(C.URL_BIOMART, params={"query": xml}, timeout=C.REQUEST_TIMEOUT)
        r.raise_for_status()
        for row in csv.reader(io.StringIO(r.text), delimiter="\t"):
            if not row or not row[0]:
                continue
            hsa, o_gene, o_name, o_type, o_perc = (row + [""] * 5)[:5]
            if not o_gene:
                continue
            # keep the best ortholog per human gene (prefer one2one, then %id)
            rec = {"ortholog_ensembl": o_gene, "ortholog_symbol": o_name or None,
                   "orthology_type": o_type or None,
                   "perc_id": float(o_perc) if o_perc else None}
            prev = out.get(hsa)
            if prev is None or _better(rec, prev):
                out[hsa] = rec
        time.sleep(0.2)
    return out


def _better(a, b):
    rank = {"ortholog_one2one": 0, "ortholog_one2many": 1, "ortholog_many2many": 2}
    ra, rb = rank.get(a["orthology_type"], 3), rank.get(b["orthology_type"], 3)
    if ra != rb:
        return ra < rb
    return (a["perc_id"] or 0) > (b["perc_id"] or 0)


def run():
    unified = json.loads((C.DATA_PROC / "unified_core.json").read_text())
    ens_ids = sorted({g["ensembl_gene"] for g in unified.values() if g.get("ensembl_gene")})

    mouse = fetch_orthologs(ens_ids, "mouse")
    fish = fetch_orthologs(ens_ids, "zebrafish")

    cons = {}
    for gid, g in unified.items():
        e = g.get("ensembl_gene")
        cons[gid] = {"mouse": mouse.get(e), "zebrafish": fish.get(e)}

    (C.DATA_PROC / "conservation.json").write_text(json.dumps(cons, indent=2))
    n_mouse = sum(1 for v in cons.values() if v["mouse"])
    n_fish = sum(1 for v in cons.values() if v["zebrafish"])
    n_either = sum(1 for v in cons.values() if v["mouse"] or v["zebrafish"])
    n = len(cons)
    print(f"[step4] mouse ortholog={n_mouse} ({100*n_mouse/n:.1f}%) "
          f"zebrafish={n_fish} ({100*n_fish/n:.1f}%) "
          f"either={n_either} ({100*n_either/n:.1f}%)")
    return cons


if __name__ == "__main__":
    run()
