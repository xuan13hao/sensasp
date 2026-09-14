"""Step 1 — Acquire the four senescence/SASP seed gene lists.

Downloads CellAge, GenAge, SenMayo and the Reactome Cellular Senescence pathway
participants, extracts gene symbols from each, and writes seed_lists_raw.json.
"""
import io
import json
import zipfile
import csv
import re

import requests

from . import config as C


def _get(url, **kw):
    r = requests.get(url, timeout=C.REQUEST_TIMEOUT, **kw)
    r.raise_for_status()
    return r


def fetch_cellage():
    """CellAge: literature-curated + predicted senescence genes (TSV inside a zip)."""
    z = zipfile.ZipFile(io.BytesIO(_get(C.URL_CELLAGE).content))
    name = next(n for n in z.namelist() if n.endswith(".tsv"))
    text = z.read(name).decode("utf-8", errors="replace")
    genes = set()
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        s = (row.get("Gene symbol") or "").strip()
        if s:
            genes.add(s)
    return sorted(genes)


def fetch_genage():
    """GenAge human: curated ageing genes (CSV inside a zip)."""
    z = zipfile.ZipFile(io.BytesIO(_get(C.URL_GENAGE).content))
    name = next(n for n in z.namelist() if n.endswith("genage_human.csv"))
    text = z.read(name).decode("utf-8-sig", errors="replace")
    genes = set()
    for row in csv.DictReader(io.StringIO(text)):
        s = (row.get("symbol") or "").strip()
        if s:
            genes.add(s)
    return sorted(genes)


def fetch_senmayo():
    """SenMayo (SAUL_SEN_MAYO, MSigDB M45803): curated senescence/SASP signature."""
    data = _get(C.URL_SENMAYO).json()
    return sorted(set(data["SAUL_SEN_MAYO"]["geneSymbols"]))


def fetch_reactome():
    """Reactome Cellular Senescence (R-HSA-2559583): gene participants.

    Parses reference entities of every physical entity in the pathway subtree.
    """
    parts = _get(C.URL_REACTOME_PARTS).json()
    genes = set()
    for pe in parts:
        for ref in pe.get("refEntities", []):
            sc = ref.get("schemaClass", "")
            disp = ref.get("displayName", "")
            ident = ref.get("identifier", "")
            if sc == "ReferenceDNASequence" and ident.startswith("ENSG"):
                genes.add(disp.split()[-1])
            elif sc == "ReferenceGeneProduct" and disp:
                genes.add(disp.split()[-1])
    # keep plausible gene symbols only
    return sorted(g for g in genes if re.match(r"^[A-Za-z0-9\-\.]+$", g) and not g.startswith("ENSG"))


def run():
    seed = {
        "CellAge": {
            "description": "HAGR CellAge — curated + predicted cellular senescence genes",
            "url": "https://genomics.senescence.info/cells/",
            "genes": fetch_cellage(),
        },
        "GenAge": {
            "description": "HAGR GenAge human — curated ageing genes",
            "url": "https://genomics.senescence.info/genes/",
            "genes": fetch_genage(),
        },
        "SenMayo": {
            "description": "SenMayo (SAUL_SEN_MAYO, MSigDB M45803) — Saul et al. 2022",
            "url": "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SAUL_SEN_MAYO",
            "genes": fetch_senmayo(),
        },
        "Reactome": {
            "description": f"Reactome Cellular Senescence pathway {C.REACTOME_PATHWAY_ID}",
            "url": f"https://reactome.org/content/detail/{C.REACTOME_PATHWAY_ID}",
            "genes": fetch_reactome(),
        },
    }
    out = C.DATA_RAW / "seed_lists_raw.json"
    out.write_text(json.dumps(seed, indent=2))
    n_total = len(set().union(*[set(v["genes"]) for v in seed.values()]))
    print(f"[step1] CellAge={len(seed['CellAge']['genes'])} "
          f"GenAge={len(seed['GenAge']['genes'])} "
          f"SenMayo={len(seed['SenMayo']['genes'])} "
          f"Reactome={len(seed['Reactome']['genes'])} | "
          f"{n_total} unique symbols -> {out}")
    return seed


if __name__ == "__main__":
    run()
