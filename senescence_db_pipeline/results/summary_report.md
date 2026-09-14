# Unified Cellular Senescence / SASP Gene Database — Build Report

## Seed sources
| Source | Genes |
|---|---|
| CellAge | 866 |
| GenAge | 307 |
| SenMayo | 124 |
| Reactome | 163 |
| **Summed** | **1460** |

## RQ1 — Unification & de-duplication
- Unique genes after unification: **1250**
- Redundant entries collapsed: **210** (14.4%)
- Genes confirmed by >=2 sources: **173**
- n_sources distribution: {'1': 1077, '2': 138, '3': 33, '4': 2}

## RQ2 — Net-new annotation coverage
| Layer | Genes | Coverage |
|---|---|---|
| Mouse ortholog | 1180 | 94.4% |
| Zebrafish ortholog | 1045 | 83.6% |
| Conservation (either) | 1198 | 95.8% |
| GTEx expression | 1224 | 97.9% |
| HPA specificity | 1224 | 97.9% |
| STRING interaction | 1162 | 93.0% |
| **All layers complete** | **1117** | **89.4%** |

None of the four seed sources provides conservation, expression or interaction
annotation, so all of the above is net-new relative to the inputs.
