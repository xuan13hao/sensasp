"""Step 8 — Compute outcome metrics (RQ1/RQ2) and write the summary report.

RQ1 = redundancy removed by unification; RQ2 = net-new annotation coverage from
the three layers. Writes build_metrics.json and summary_report.md.
"""
import json

from . import config as C


def run():
    seed = json.loads((C.DATA_RAW / "seed_lists_raw.json").read_text())
    db = json.loads((C.RESULTS / "senescence_sasp_database.json").read_text())
    genes = db["genes"]
    n = len(genes)

    summed = sum(len(seed[s]["genes"]) for s in C.SOURCES)
    redundant = summed - n

    # RQ1: cross-source confirmation
    by_ns = {}
    for g in genes.values():
        by_ns[g["n_sources"]] = by_ns.get(g["n_sources"], 0) + 1
    multi = sum(v for k, v in by_ns.items() if k >= 2)

    # RQ2: layer coverage
    n_mouse = sum(1 for g in genes.values() if (g["conservation"] or {}).get("mouse"))
    n_fish = sum(1 for g in genes.values() if (g["conservation"] or {}).get("zebrafish"))
    n_cons = sum(1 for g in genes.values()
                 if (g["conservation"] or {}).get("mouse") or (g["conservation"] or {}).get("zebrafish"))
    n_gtex = sum(1 for g in genes.values() if (g["expression"] or {}).get("gtex_median_tpm"))
    n_hpa = sum(1 for g in genes.values() if (g["expression"] or {}).get("hpa"))
    n_int = sum(1 for g in genes.values() if (g["interaction"] or {}).get("string_degree", 0) > 0)
    n_all = sum(1 for g in genes.values()
                if ((g["conservation"] or {}).get("mouse") or (g["conservation"] or {}).get("zebrafish"))
                and (g["expression"] or {}).get("gtex_median_tpm")
                and (g["interaction"] or {}).get("string_degree", 0) > 0)

    def pct(x):
        return round(100 * x / n, 1)

    metrics = {
        "sources": {s: len(seed[s]["genes"]) for s in C.SOURCES},
        "summed_entries": summed,
        "unique_genes": n,
        "redundant_collapsed": redundant,
        "redundancy_pct": round(100 * redundant / summed, 1),
        "n_sources_distribution": {str(k): by_ns[k] for k in sorted(by_ns)},
        "multi_source_genes": multi,
        "coverage": {
            "mouse_ortholog": {"n": n_mouse, "pct": pct(n_mouse)},
            "zebrafish_ortholog": {"n": n_fish, "pct": pct(n_fish)},
            "conservation_any": {"n": n_cons, "pct": pct(n_cons)},
            "gtex_expression": {"n": n_gtex, "pct": pct(n_gtex)},
            "hpa_specificity": {"n": n_hpa, "pct": pct(n_hpa)},
            "string_interaction": {"n": n_int, "pct": pct(n_int)},
            "all_layers_complete": {"n": n_all, "pct": pct(n_all)},
        },
    }
    (C.RESULTS / "build_metrics.json").write_text(json.dumps(metrics, indent=2))

    md = f"""# Unified Cellular Senescence / SASP Gene Database — Build Report

## Seed sources
| Source | Genes |
|---|---|
""" + "\n".join(f"| {s} | {len(seed[s]['genes'])} |" for s in C.SOURCES) + f"""
| **Summed** | **{summed}** |

## RQ1 — Unification & de-duplication
- Unique genes after unification: **{n}**
- Redundant entries collapsed: **{redundant}** ({metrics['redundancy_pct']}%)
- Genes confirmed by >=2 sources: **{multi}**
- n_sources distribution: {metrics['n_sources_distribution']}

## RQ2 — Net-new annotation coverage
| Layer | Genes | Coverage |
|---|---|---|
| Mouse ortholog | {n_mouse} | {pct(n_mouse)}% |
| Zebrafish ortholog | {n_fish} | {pct(n_fish)}% |
| Conservation (either) | {n_cons} | {pct(n_cons)}% |
| GTEx expression | {n_gtex} | {pct(n_gtex)}% |
| HPA specificity | {n_hpa} | {pct(n_hpa)}% |
| STRING interaction | {n_int} | {pct(n_int)}% |
| **All layers complete** | **{n_all}** | **{pct(n_all)}%** |

None of the four seed sources provides conservation, expression or interaction
annotation, so all of the above is net-new relative to the inputs.
"""
    (C.RESULTS / "summary_report.md").write_text(md)
    print(f"[step8] unique={n} redundant={redundant} ({metrics['redundancy_pct']}%) "
          f"all_layers={n_all} ({pct(n_all)}%)")
    return metrics


if __name__ == "__main__":
    run()
