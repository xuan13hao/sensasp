#!/usr/bin/env python3
"""
SenSASP MCP Server
Supports both stdio (local) and SSE/HTTP (remote) transports.

Usage:
  sensasp-mcp                        # stdio (Claude Desktop default)
  sensasp-mcp --transport sse        # HTTP server on port 8080
  sensasp-mcp --transport sse --port 3000
"""

import argparse
import json
from importlib.resources import files

from mcp.server.fastmcp import FastMCP

# ── Load bundled data ──────────────────────────────────────────────────────
_pkg = files("sensasp_mcp")

_FLAT: list[dict]  = json.loads((_pkg / "data" / "genes.json").read_text())
_FULL: list[dict]  = json.loads((_pkg / "data" / "genes_full.json").read_text())
_META: dict        = json.loads((_pkg / "data" / "metrics.json").read_text())
_SOURCES: list     = json.loads((_pkg / "data" / "sources.json").read_text())

_FULL_INDEX: dict[str, dict] = {
    g["identifiers"]["hgnc_symbol"].upper(): g for g in _FULL
}

# Secondary index of original source (input) symbols so queries using the symbol
# as it appeared in a seed resource (e.g. 'ATP5O' → canonical 'ATP5PO') still
# resolve. Canonical symbols take precedence; ambiguous source symbols (mapping
# to more than one canonical record) are dropped to avoid wrong matches.
def _build_source_index(records: list[dict]) -> dict[str, dict]:
    index: dict[str, dict] = {}
    ambiguous: set[str] = set()
    for g in records:
        for alias in g["identifiers"].get("input_symbols") or []:
            a = alias.upper()
            if a in _FULL_INDEX:            # canonical symbol wins
                continue
            if a in index and index[a] is not g:
                ambiguous.add(a)
            else:
                index[a] = g
    for a in ambiguous:
        index.pop(a, None)
    return index

_SOURCE_INDEX: dict[str, dict] = _build_source_index(_FULL)

# canonical symbol → original source (input) symbols that differ from it, so
# free-text search can also match a gene by the symbol it carried in a seed
# resource (e.g. searching 'ATP5O' surfaces canonical 'ATP5PO').
_ALIASES: dict[str, list[str]] = {}
for _g in _FULL:
    _canon = _g["identifiers"]["hgnc_symbol"].upper()
    _extra = [s for s in (_g["identifiers"].get("input_symbols") or []) if s.upper() != _canon]
    if _extra:
        _ALIASES[_canon] = _extra

# ── MCP server ─────────────────────────────────────────────────────────────
mcp = FastMCP(
    "SenSASP",
    instructions=(
        "SenSASP is a unified, multi-layer gene database for cellular senescence "
        "and the senescence-associated secretory phenotype (SASP). It integrates "
        "1,250 unique human genes from CellAge, GenAge, SenMayo, and the Reactome "
        "Cellular Senescence pathway, annotated with cross-species conservation "
        "(mouse & zebrafish), GTEx tissue expression (54 tissues), HPA single-cell "
        "specificity, and STRING protein-protein interactions (score >= 0.700). "
        "Cite: Xuan et al. 2026, bioRxiv doi:10.64898/2026.08.26.747427"
    ),
)


@mcp.tool()
def search_genes(
    query: str = "",
    source: str = "",
    multi_source_only: bool = False,
    conserved_mouse: bool = False,
    conserved_zebrafish: bool = False,
    hpa_tissue_specificity: str = "",
    hpa_cell_specificity: str = "",
    min_string_degree: int = 0,
    limit: int = 20,
) -> list[dict]:
    """
    Search the SenSASP gene database with optional filters.

    Args:
        query: Free-text against symbol, original source (input) symbols,
            gene name, and STRING partners.
        source: Filter by seed source — CellAge, GenAge, SenMayo, or Reactome.
        multi_source_only: Only return genes confirmed by >=2 sources.
        conserved_mouse: Only return genes with a mouse ortholog.
        conserved_zebrafish: Only return genes with a zebrafish ortholog.
        hpa_tissue_specificity: HPA RNA *tissue* specificity category, e.g.
            'Tissue enhanced', 'Tissue enriched', 'Low tissue specificity'.
        hpa_cell_specificity: HPA RNA *single-cell* specificity category, e.g.
            'Cell type enhanced', 'Cell type enriched', 'Low cell type specificity'.
        min_string_degree: Minimum STRING intra-set interaction degree.
        limit: Max results (default 20, max 200).
    """
    limit = min(limit, 200)
    q = query.lower()
    exact, others = [], []
    for g in _FLAT:
        aliases = _ALIASES.get(g["symbol"].upper(), [])
        if q and q not in f"{g['symbol']} {g['gene_name']} {g['string_top_partners']} {' '.join(aliases)}".lower():
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
        if hpa_cell_specificity and g["hpa_cell_spec"] != hpa_cell_specificity:
            continue
        if min_string_degree and (g["string_degree"] or 0) < min_string_degree:
            continue
        row = {
            "symbol":              g["symbol"],
            "input_symbols":       aliases,
            "gene_name":           g["gene_name"],
            "sources":             g["sources"],
            "n_sources":           g["n_sources"],
            "ensembl":             g["ensembl"],
            "conserved_mouse":     g["conserved_mouse"],
            "conserved_zebrafish": g["conserved_zebrafish"],
            "gtex_top_tissue":     g["gtex_top_tissue"],
            "gtex_max_tpm":        g["gtex_max_tpm"],
            "hpa_tissue_spec":     g["hpa_tissue_spec"],
            "hpa_cell_spec":       g["hpa_cell_spec"],
            "string_degree":       g["string_degree"],
            "string_top_partners": g["string_top_partners"],
        }
        # Exact symbol match (canonical or source symbol) floats to the top
        if q and (g["symbol"].lower() == q or any(a.lower() == q for a in aliases)):
            exact.append(row)
        else:
            others.append(row)
    return (exact + others)[:limit]


@mcp.tool()
def get_gene(symbol: str) -> dict:
    """
    Full annotation record for a single gene.

    Args:
        symbol: HGNC gene symbol, e.g. 'TP53', 'IL6', 'CDKN2A'. Case-insensitive.

    Returns full provenance, conservation, GTEx top tissues, HPA specificity,
    and STRING partners with scores. Original source symbols (as they appeared
    in the seed resources) also resolve to the canonical record when unambiguous.
    """
    key = symbol.upper()
    record = _FULL_INDEX.get(key) or _SOURCE_INDEX.get(key)
    if record is None:
        return {"error": f"Gene '{symbol}' not found in SenSASP (1,250 genes)."}
    return record


@mcp.tool()
def get_database_stats() -> dict:
    """
    Overall statistics and annotation coverage for the SenSASP database.
    """
    return {
        "name":                   "SenSASP",
        "build_date":             _META["build_date"],
        "n_genes":                _META["n_genes"],
        "sources":                _META["sources"],
        "n_sources_distribution": _META["n_sources_distribution"],
        "multi_source_genes":     _META["multi_source_genes"],
        "coverage":               _META["coverage"],
        "web_url":                "https://xuan13hao.github.io/sensasp",
        "doi":                    "10.64898/2026.08.26.747427",
    }


@mcp.tool()
def list_sources() -> list[dict]:
    """
    The four seed databases integrated in SenSASP with gene counts and URLs.
    """
    return _SOURCES


@mcp.tool()
def get_top_hub_genes(n: int = 20) -> list[dict]:
    """
    Top-N genes by STRING intra-set interaction degree (hub genes).

    Args:
        n: Number to return (default 20, max 100).
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
    All genes listed by 2 or more seed sources, sorted by source overlap (n_sources desc).

    n_sources counts how many seed resources list each gene, not experimental
    strength. The seed resources use different inclusion criteria, so higher overlap
    means broader agreement across databases rather than stronger confidence on its
    own. These 173 genes are the most-shared across databases; consult the original
    studies (see list_sources) when assessing experimental evidence for a gene.
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
def main():
    parser = argparse.ArgumentParser(description="SenSASP MCP Server")
    parser.add_argument(
        "--transport",
        choices=["stdio", "sse"],
        default="stdio",
        help="Transport type: stdio (local/Claude Desktop) or sse (remote HTTP)",
    )
    parser.add_argument(
        "--port", type=int, default=8080,
        help="Port for SSE transport (default: 8080)",
    )
    parser.add_argument(
        "--host", default="0.0.0.0",
        help="Host for SSE transport (default: 0.0.0.0)",
    )
    args = parser.parse_args()

    if args.transport == "sse":
        print(f"SenSASP MCP server running on http://{args.host}:{args.port}/sse")
        mcp.run(transport="sse", port=args.port, host=args.host)
    else:
        mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
