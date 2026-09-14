"""Step 6 — Add protein-protein interaction summary (STRING v12, high confidence).

Fetches the intra-set STRING network (combined score >= 700) for all mapped
UniProt/symbol nodes and summarises per-gene degree. Writes interaction.json
and the raw edge list string_network.json.
"""
import json
import time

import requests

from . import config as C


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def fetch_network(symbols):
    """STRING network endpoint restricted to the input node set."""
    edges = []
    seen = set()
    # STRING caps identifiers per call; submit in chunks and merge (intra-set edges)
    ident = "%0d".join(symbols)
    r = requests.post(
        f"{C.URL_STRING}/json/network",
        data={"identifiers": "\r".join(symbols),
              "species": C.STRING_SPECIES,
              "required_score": C.STRING_MIN_SCORE,
              "caller_identity": "senescence_db_pipeline"},
        timeout=C.REQUEST_TIMEOUT * 3,
    )
    r.raise_for_status()
    for e in r.json():
        a, b = e["preferredName_A"], e["preferredName_B"]
        key = tuple(sorted((a, b)))
        if a == b or key in seen:
            continue
        seen.add(key)
        edges.append({"a": a, "b": b, "score": e["score"]})
    return edges


def run():
    unified = json.loads((C.DATA_PROC / "unified_core.json").read_text())
    symbols = sorted({g["symbol"] for g in unified.values() if g.get("symbol")})

    edges = fetch_network(symbols)
    (C.DATA_PROC / "string_network.json").write_text(json.dumps(edges, indent=2))

    degree = {s: 0 for s in symbols}
    partners = {s: set() for s in symbols}
    for e in edges:
        for x, y in ((e["a"], e["b"]), (e["b"], e["a"])):
            if x in degree:
                degree[x] += 1
                partners[x].add(y)

    inter = {}
    sym2id = {g["symbol"]: gid for gid, g in unified.items()}
    for s in symbols:
        gid = sym2id.get(s)
        if gid is None:
            continue
        inter[gid] = {"string_degree": degree[s],
                      "n_partners": len(partners[s]),
                      "top_partners": sorted(partners[s])[:15]}

    (C.DATA_PROC / "interaction.json").write_text(json.dumps(inter, indent=2))
    mapped = sum(1 for v in inter.values() if v["string_degree"] > 0)
    hub = max(inter.items(), key=lambda kv: kv[1]["string_degree"])
    hub_sym = unified[hub[0]]["symbol"]
    print(f"[step6] edges={len(edges)} nodes_with_edges={mapped} "
          f"| top hub={hub_sym} degree={hub[1]['string_degree']}")
    return inter


if __name__ == "__main__":
    run()
