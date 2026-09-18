#!/usr/bin/env python3
"""
SenSASP MCP Server
Exposes the SenSASP unified senescence / SASP gene database as MCP tools.
"""

import json
from pathlib import Path
from mcp.server.fastmcp import FastMCP

# ── Load data ──────────────────────────────────────────────────────────────
_ROOT = Path(__file__).parent
_DB_PATH   = _ROOT / "senescence_sasp_database.json"
_FLAT_PATH = _ROOT / "docs" / "data" / "genes.json"
_METRICS   = _ROOT / "senescence_db_pipeline" / "results" / "build_metrics.json"

with open(_DB_PATH) as f:
    _DB = json.load(f)

with open(_FLAT_PATH) as f:
    _FLAT: list[dict] = json.load(f)

with open(_METRICS) as f:
    _METRICS_DATA = json.load(f)

# Build symbol → full record index for O(1) lookup
_FULL_INDEX: dict[str, dict] = {
    g["identifiers"]["hgnc_symbol"].upper(): g
    for g in _DB["genes"]
}

# ── Server ─────────────────────────────────────────────────────────────────
mcp = FastMCP(
    "SenSASP",
    instructions=(
        "SenSASP is a unified, multi-layer gene database for cellular senescence "
        "and the senescence-associated secretory phenotype (SASP). It integrates "
        "1,250 unique human genes from CellAge, GenAge, SenMayo, and the Reactome "
        "Cellular Senescence pathway, annotated with cross-species conservation, "
        "GTEx tissue expression, HPA single-cell specificity, and STRING "
        "protein–protein interactions. Use search_genes to find genes matching "
        "criteria, get_gene for full details on a single gene, and "
        "get_database_stats for an overview of the database."
    ),
)


# ── Tools ──────────────────────────────────────────────────────────────────

@mcp.tool()
def search_genes(
    query: str = "",
    source: str = "",
    multi_source_only: bool = False,
    conserved_mouse: bool = False,
    conserved_zebrafish: bool = False,
    hpa_tissue_specificity: str = "",
    min_string_degree: int = 0,
    limit: int = 20,
) -> list[dict]:
    """
    Search the SenSASP gene database with optional filters.

    Args:
        query: Free-text search against gene symbol, gene name, and STRING
               interaction partners (case-insensitive, partial match).
        source: Filter by seed source. One of: CellAge, GenAge, SenMayo, Reactome.
                Leave empty to include all sources.
        multi_source_only: If True, return only genes confirmed by 2 or more sources.
        conserved_mouse: If True, return only genes with a mouse ortholog.
        conserved_zebrafish: If True, return only genes with a zebrafish ortholog.
        hpa_tissue_specificity: Filter by HPA RNA tissue specificity category.
                                Examples: 'Tissue enhanced', 'Tissue specific',
                                'Low cell type specificity', 'Cell type enhanced'.
        min_string_degree: Minimum STRING intra-set interaction degree.
        limit: Maximum number of results to return (default 20, max 200).

    Returns:
        List of matching gene records with key fields.
    """
    limit = min(limit, 200)
    q = query.lower()

    results = []
    for g in _FLAT:
        if q and q not in f"{g['symbol']} {g['gene_name']} {g['string_top_partners']}".lower():
            continue
        if source and source not in g["sources"]:
            continue
        if multi_source_only and g["n_sources"] < 2:
            continue
        if conserved_mouse and not g["conserved_mouse"]:
            continue
        if conserved_zebrafish and not g["conserved_zebrafish"]:
            continue
        if hpa_tissue_specificity and g["hpa_tissue_spec"] != hpa_tissue_specificity:
            continue
        if min_string_degree and (g["string_degree"] or 0) < min_string_degree:
            continue

        results.append({
            "symbol":          g["symbol"],
            "gene_name":       g["gene_name"],
            "sources":         g["sources"],
            "n_sources":       g["n_sources"],
            "ensembl":         g["ensembl"],
            "conserved_mouse": g["conserved_mouse"],
            "conserved_zebrafish": g["conserved_zebrafish"],
            "gtex_top_tissue": g["gtex_top_tissue"],
            "gtex_max_tpm":    g["gtex_max_tpm"],
            "hpa_tissue_spec": g["hpa_tissue_spec"],
            "string_degree":   g["string_degree"],
            "string_top_partners": g["string_top_partners"],
        })
        if len(results) >= limit:
            break

    return results


@mcp.tool()
def get_gene(symbol: str) -> dict:
    """
    Retrieve the full annotation record for a single gene.

    Args:
        symbol: HGNC gene symbol (e.g. 'TP53', 'IL6', 'CDKN2A'). Case-insensitive.

    Returns:
        Full gene record including identifiers, provenance, conservation,
        GTEx expression (all top tissues), HPA specificity, and STRING
        interaction partners with scores. Returns an error dict if not found.
    """
    key = symbol.upper()
    record = _FULL_INDEX.get(key)
    if record is None:
        return {"error": f"Gene '{symbol}' not found in SenSASP database."}
    return record


@mcp.tool()
def get_database_stats() -> dict:
    """
    Return overall statistics and coverage metrics for the SenSASP database.

    Returns:
        Summary including total gene count, per-source counts, multi-source
        confirmation numbers, and annotation coverage percentages.
    """
    return {
        "name":        _DB["metadata"]["name"],
        "build_date":  _DB["metadata"]["build_date"],
        "n_genes":     _DB["metadata"]["n_genes"],
        "sources":     _METRICS_DATA["sources"],
        "n_sources_distribution": _METRICS_DATA["n_sources_distribution"],
        "multi_source_genes":     _METRICS_DATA["multi_source_genes"],
        "coverage":    _METRICS_DATA["coverage"],
        "web_url":     "https://xuan13hao.github.io/sensasp",
        "doi":         "10.64898/2026.08.26.747427",
    }


@mcp.tool()
def list_sources() -> list[dict]:
    """
    List the four seed databases integrated in SenSASP with their metadata.

    Returns:
        Each source's name, gene count, description, and URL.
    """
    return [
        {
            "name":        name,
            "n_genes":     info["n_input"],
            "description": info["description"],
            "url":         info["url"],
        }
        for name, info in _DB["metadata"]["seed_sources"].items()
    ]


@mcp.tool()
def get_top_hub_genes(n: int = 20) -> list[dict]:
    """
    Return the top-N genes by STRING intra-set interaction degree (hub genes).

    Args:
        n: Number of genes to return (default 20, max 100).

    Returns:
        Genes sorted descending by STRING degree, with symbol, sources,
        degree, and top interaction partners.
    """
    n = min(n, 100)
    ranked = sorted(
        [g for g in _FLAT if (g["string_degree"] or 0) > 0],
        key=lambda g: g["string_degree"],
        reverse=True,
    )
    return [
        {
            "symbol":              g["symbol"],
            "gene_name":           g["gene_name"],
            "sources":             g["sources"],
            "n_sources":           g["n_sources"],
            "string_degree":       g["string_degree"],
            "string_top_partners": g["string_top_partners"],
        }
        for g in ranked[:n]
    ]


@mcp.tool()
def get_multi_source_genes() -> list[dict]:
    """
    Return all genes confirmed by 2 or more seed sources, sorted by n_sources descending.

    Returns:
        Genes corroborated by multiple databases, indicating higher confidence
        senescence relevance. Includes which sources confirm each gene.
    """
    multi = sorted(
        [g for g in _FLAT if g["n_sources"] >= 2],
        key=lambda g: g["n_sources"],
        reverse=True,
    )
    return [
        {
            "symbol":    g["symbol"],
            "gene_name": g["gene_name"],
            "n_sources": g["n_sources"],
            "sources":   g["sources"],
            "ensembl":   g["ensembl"],
        }
        for g in multi
    ]


# ── Entry point ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    mcp.run()
