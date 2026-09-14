# A Unified, Multi-Layer Reference Database of Cellular Senescence and SASP Genes

## Abstract

Research on cellular senescence and the senescence-associated secretory phenotype
(SASP) draws on several independently curated gene resources that use different
scopes, identifier conventions, and update cycles, making cross-resource
integration error-prone. We unified four widely used seed resources — CellAge,
GenAge (human), the SenMayo signature, and the Reactome Cellular Senescence
pathway — onto a single canonical identifier space and enriched every gene with
three annotation layers absent from all four inputs: cross-species conservation
(mouse and zebrafish orthologs), tissue and cell-type expression (GTEx v8 and the
Human Protein Atlas), and high-confidence protein–protein interactions (STRING
v12). Unification collapsed 1,460 summed source entries into 1,250 unique genes
(210 redundant entries removed, 14.4%), while preserving full source provenance:
173 genes are corroborated by two or more resources and two (IL6, JUN) appear in
all four. The annotation layers reach 95.8% of genes for conservation, 97.8% for
expression, and 92.9% for interactions, with 89.3% of genes annotated across all
three layers. The result is a single machine-readable database with harmonized
identifiers and net-new functional context, released together with a fully
scripted, one-command-reproducible pipeline.

## 1. Introduction

Cellular senescence — a durable arrest of proliferation accompanied by a
pro-inflammatory secretome (the SASP) — is central to ageing, tissue remodeling,
tumor suppression, and age-related disease. As the field has grown, the community
has produced several curated gene resources, each reflecting a different
definition and evidence base. CellAge catalogs genes whose manipulation drives or
prevents senescence; GenAge collects genes implicated in ageing more broadly;
SenMayo is an experimentally derived SASP/senescence transcriptomic signature;
and the Reactome Cellular Senescence pathway encodes a curated mechanistic model.

These resources are complementary but disconnected. They differ in scope (a
senescence-specific screen versus a general ageing catalog versus a secretory
signature), in identifier conventions (gene symbols that drift over time versus
pathway participant records versus signature membership), and in update cadence.
A researcher who wants "the senescence gene set" must therefore reconcile four
files by hand, resolving outdated symbols and duplicate entries, with no shared
key to join them. Crucially, none of the four provides the functional context
that most downstream analyses require: whether a gene is conserved in common
model organisms, where it is expressed, and how it is wired into the interaction
network.

We address both gaps. First, we unify the four resources onto a single canonical
identifier (Ensembl gene ID), removing redundancy while recording exactly which
sources support each gene (research question RQ1). Second, we add three
annotation layers — conservation, expression, and interactions — that are absent
from every input, quantifying the net-new coverage they provide (RQ2). The
deliverable is a single JSON/CSV database plus a scripted pipeline that
regenerates it end to end from public APIs.

## 2. Methods

### 2.1 Seed acquisition

Four seed lists were downloaded programmatically from their primary
distributors. **CellAge** and **GenAge (human)** were obtained from the Human
Ageing Genomic Resources (HAGR) as versioned archives and parsed for the curated
gene-symbol columns. **SenMayo** was retrieved as the `SAUL_SEN_MAYO` gene set
(MSigDB accession M45803; Saul et al., 2022) via the MSigDB JSON download
endpoint. The **Reactome Cellular Senescence** pathway (R-HSA-2559583) was
retrieved through the Reactome ContentService participants endpoint; gene-level
participants were extracted from the reference entities of every physical entity
in the pathway subtree (`ReferenceDNASequence` records carrying Ensembl gene IDs
and `ReferenceGeneProduct` records carrying gene symbols).

### 2.2 Identifier harmonization

Every unique input symbol was mapped to a canonical identifier set — current HGNC
symbol, Ensembl gene ID, UniProt accession, Entrez ID, and HGNC ID — using
MyGene.info with alias-aware scopes (`symbol,alias`) restricted to human. When a
query returned multiple hits, we preferred human records with an exact
current-symbol match and an assigned Ensembl gene, breaking remaining ties by
MyGene relevance score. Symbols that MyGene resolved to a different current
symbol were flagged as updated.

### 2.3 Unification (RQ1)

Harmonized genes were collapsed onto the Ensembl gene ID as the canonical key
(falling back to current symbol, then input symbol, for the small number of
non-coding entries lacking an Ensembl mapping). For each unified gene we recorded
the set of contributing sources and the set of original input symbols, yielding a
per-gene source-provenance vector used for all overlap statistics.

### 2.4 Conservation layer

Mouse (*Mus musculus*) and zebrafish (*Danio rerio*) orthologs were retrieved
from Ensembl BioMart (martservice REST/XML) in chunks of 200 human Ensembl gene
IDs, requesting ortholog gene ID, ortholog symbol, orthology type, and percent
identity. Where multiple orthologs existed, the best was selected by orthology
type (one2one > one2many > many2many) then percent identity.

### 2.5 Expression layer

Human Ensembl gene IDs were resolved to versioned GENCODE IDs through the GTEx
reference endpoint, then median TPM values were pulled across all 54 GTEx v8
tissues from the median-gene-expression endpoint. Independently, the Human
Protein Atlas `search_download` API supplied RNA tissue-specificity class, RNA
tissue distribution, single-cell-type specificity, and subcellular location per
gene.

### 2.6 Interaction layer

High-confidence protein–protein interactions were obtained from STRING v12 for
the full set of mapped gene symbols (*Homo sapiens*, taxon 9606) using a combined
score threshold of 700. The returned network was restricted to intra-set edges,
and per-gene degree and interaction partners were computed from the resulting
graph.

### 2.7 Assembly, metrics, and figures

The unified core and the three annotation layers were merged into a single
nested JSON database keyed by canonical ID, with a flattened CSV view for
spreadsheet use. Outcome metrics for RQ1 (redundancy removed, cross-source
confirmation) and RQ2 (per-layer coverage) were computed from the assembled
database. All figures were generated from the pipeline outputs in a single
Jupyter notebook.

### 2.8 Reproducibility

The pipeline is implemented as eight ordered Python modules (`step1`–`step8`)
driven by `run_pipeline.py`, with a `main.sh` wrapper that installs dependencies,
runs the data pipeline, and executes the figure notebook end to end. All data are
fetched from public APIs at run time; no manual downloads are required.

## 3. Data

The database is distributed as `senescence_sasp_database.json` (nested, one record
per gene) and `senescence_sasp_database.csv` (flat view). Each gene record
contains: the canonical Ensembl ID; current HGNC symbol and gene name; an
identifier block (Ensembl, UniProt, Entrez, HGNC ID, and all original input
symbols); the source-provenance list and source count; and the three annotation
blocks. The **conservation** block gives mouse and zebrafish ortholog ID, symbol,
orthology type, and percent identity. The **expression** block gives the GENCODE
ID, per-tissue GTEx median TPM across 54 tissues, and the HPA specificity and
localization fields. The **interaction** block gives STRING degree, partner
count, and top partners.

Provenance is preserved end to end: the raw seed lists
(`seed_lists_raw.json`), the full identifier map (`id_map.csv`), and the
source-overlap matrix (`source_overlap.csv`) are all retained as intermediate
artifacts, so any gene can be traced back to the exact symbols and sources it was
built from.

## 4. Results

### 4.1 Unification removes substantial redundancy while preserving provenance (RQ1)

The four seed lists contributed 866 (CellAge), 307 (GenAge), 124 (SenMayo), and
163 (Reactome) genes — 1,460 entries summed across sources, spanning 1,252 unique
input symbols. Identifier harmonization mapped all 1,252 symbols to Ensembl gene
IDs and 1,249 to UniProt accessions (the three unmapped-to-UniProt cases are
non-coding genes such as TERC), and updated 5 outdated symbols to current HGNC
nomenclature (e.g., ARNTL→BMAL1, CTGF→CCN2, H2AFX→H2AX).

Collapsing onto the canonical key produced **1,250 unique genes**, removing **210
redundant entries (14.4%)** relative to the summed input. Source provenance shows
that the great majority of genes are resource-specific: 1,077 genes come from a
single source, while **173 are corroborated by two or more** resources — 138 by
exactly two, 33 by three, and **2 (IL6 and JUN) by all four**. This distribution
quantifies both the fragmentation of the current landscape and the small,
high-confidence core the resources agree on (Figure 1).

### 4.2 Three annotation layers add functional context absent from every source (RQ2)

None of the four seed resources provides conservation, expression, or interaction
annotation, so every value in these layers is net-new. Coverage is high across
all three (Figure 4):

- **Conservation:** 94.4% of genes have a mouse ortholog and 83.6% a zebrafish
  ortholog; **95.8% have an ortholog in at least one** of the two model organisms.
- **Expression:** GTEx v8 median expression across 54 tissues is available for
  97.8% of genes, and HPA specificity annotation for 97.8%.
- **Interactions:** 92.9% of genes participate in at least one high-confidence
  (score ≥ 700) STRING interaction within the set, forming a network of 15,611
  edges. The network is strongly hub-organized, with TP53 the top hub (degree
  368), consistent with its central role in senescence (Figure 3).

**89.3% of genes are annotated across all three layers simultaneously**, giving a
large, uniformly characterized core suitable for integrative analysis.

### 4.3 Expression structure of the high-confidence core

Restricting to genes supported by three or more sources and clustering their GTEx
profiles reveals coherent tissue-expression modules (Figure 2), including a set
of broadly and highly expressed regulators alongside genes with more
tissue-restricted patterns — structure that is only visible once expression is
joined to the unified gene set.

### 4.4 Application Case 1: Database-guided multi-evidence gene prioritization

A central motivation for integrating multiple evidence streams is to enable
systematic prioritization of candidate genes for experimental follow-up — a task
that is impossible from any single source alone. To illustrate this utility, we
defined a composite priority score that combines three orthogonal evidence
dimensions available in the database: (i) **source confirmation** (n_sources × 25,
maximum 100), which reflects the degree of cross-resource consensus; (ii)
**network centrality** (STRING degree normalized to the set maximum, scaled to 50),
capturing functional connectivity and likely regulatory importance; and (iii)
**evolutionary conservation** (mean mouse/*Danio rerio* sequence identity scaled by
0.25), providing independent support for biological essentiality. The three
components are intentionally heterogeneous — they capture evidence from
literature-derived curation, protein interaction topology, and comparative
genomics, respectively — and their linear combination produces a ranked candidate
list that no individual evidence layer could supply on its own.

Applied to all 1,250 unified genes, the top 20 (Figure 5) are exclusively genes
confirmed by three or four sources (score range 100–142), reflecting the strong
upward contribution of cross-source confirmation. **JUN** ranks first (score 142),
combining all-four-source confirmation (100 points) with substantial network
connectivity (21 points, degree 155) and good conservation (21 points, mouse
96.4% identity). **TP53** ranks second (141), where its exceptional network
centrality (50 points; degree 368, the highest in the set) compensates for its
lower cross-species sequence conservation (16 points), consistent with the
well-established role of p53 as the master regulator of stress responses including
oncogene-induced senescence (Vousden & Lane, 2007; Rufini et al., 2013). **IL6**
ranks third (132), reflecting its status as the canonical SASP cytokine confirmed
by all four resources (Coppe et al., 2008; Acosta et al., 2013). The remaining
top-20 entries — including CTNNB1, STAT3, EGFR, CDKN2A, ATM, and MAPK14 — are
each confirmed by three sources and represent established nodes in senescence
signaling (Campisi & d'Adda di Fagagna, 2007; Herranz & Gil, 2018).

Crucially, this ranking is only possible because the database co-locates source
membership, network degree, and conservation for every gene. A researcher using
any single seed resource cannot recover the cross-source confirmation dimension;
a researcher using STRING alone has no access to ortholog identity; and neither
provides the unified identifier space needed to join the two. The priority score
is illustrative rather than prescriptive — users may weight the components
differently, add further evidence layers (e.g., GWAS loci, single-cell
co-expression), or restrict to a particular source subset — but the underlying
database makes all such analyses a query rather than a data-engineering project.

### 4.5 Application Case 2: Aging tissue expression landscape of top-priority
senescence genes

Beyond gene prioritization, the database enables direct inspection of where
high-confidence senescence genes are expressed across tissues relevant to the
biology of ageing. Ageing is a systemic process with well-documented
tissue-heterogeneous manifestations: senescence burden accumulates differentially
in aorta, adipose tissue, liver, brain, and skeletal muscle, among other sites,
and this heterogeneity is thought to underlie the tissue-specific pathologies of
ageing (López-Otín et al., 2013; Childs et al., 2017; Birch & Gil, 2020).
Charting the expression landscape of the most evidence-supported senescence genes
across these tissues is therefore directly relevant to understanding which
processes are active in which ageing contexts, and to selecting appropriate model
systems for experimental work.

We selected 15 aging-relevant tissues from GTEx v8 — covering brain (cortex and
hippocampus), cardiovascular (heart left ventricle, aorta, coronary artery),
metabolic (liver, adipose subcutaneous and visceral), musculoskeletal (skeletal
muscle), integumentary (sun-exposed skin), renal (kidney cortex), gastrointestinal
(colon, small intestine), and haematopoietic (whole blood) compartments — and
plotted the log2(median TPM + 1) profiles for the top 15 priority-scored genes,
with rows and columns independently hierarchically clustered (Figure 6).

The resulting heatmap reveals three functionally coherent expression modules.
First, a group comprising **CTNNB1**, **STAT3**, and **RELA** shows broad and
relatively uniform expression across all tissues examined, consistent with their
roles as constitutive transcriptional regulators whose senescence-relevant activity
is controlled post-translationally rather than at the transcriptional level
(Clevers & Nusse, 2012; Yu et al., 2009). Second, **JUN** and **FOS**, the
canonical AP-1 components, display high expression particularly in adipose and
vascular tissues — sites where senescence is particularly prominent in
aged individuals and where AP-1 has been directly linked to the SASP transcriptional
program (Herranz & Gil, 2018; Lopes-Paciencia et al., 2019). Third, **IL6**
shows a distinct pattern with notably elevated expression in whole blood relative
to solid tissues, reflecting the secreted nature of IL-6 and its principal
production by immune and stromal cells — a tissue-specific expression signature
that would not be apparent from network topology or source confirmation alone
(Coppe et al., 2008). Genes with generally low median TPM in bulk tissue (ATM,
CCNA2, E2F1) reflect their cell-cycle-coupled regulation: they are expressed in
proliferating cells that are a minority in most adult tissues and become
essentially absent in post-mitotic compartments such as mature skeletal muscle
and brain neurons.

Together, the two application cases demonstrate that the value of the database
extends beyond curation and de-duplication. By co-localizing source provenance,
network topology, evolutionary conservation, and tissue expression in a unified,
queryable structure, it enables analyses — systematic prioritization, integrative
scoring, and multi-tissue expression contextualization — that would otherwise
require the bespoke integration of four separate databases, an annotation pipeline,
and the manual resolution of identifier conflicts.

## 5. Discussion

We built a single reference database that resolves the two practical problems with
the current senescence gene landscape: fragmentation across four
differently-scoped resources with drifting identifiers, and the absence of
functional context needed for downstream work. Unification onto a canonical
identifier removed 14.4% redundancy while retaining full source provenance, and
the three annotation layers added conservation, expression, and interaction
context to roughly nine of every ten genes.

The source-overlap structure is informative in itself. That only two genes (IL6,
JUN) are shared by all four resources, and only 173 by two or more, reflects
genuinely different operational definitions — a senescence-driver screen, a
general ageing catalog, a secretory signature, and a curated mechanistic pathway
capture different but overlapping biology. Preserving the per-gene source vector
lets users tune stringency to their question: a high-confidence multi-source core
for hypothesis generation, or the full union for discovery.

Several limitations follow from the design choices. We retained GenAge in full
rather than applying a strict senescence-only evidence filter, which would have
reduced it to ~25 genes; we judged that overly aggressive and instead expose the
source vector so users can filter themselves. The canonical key is the human
Ensembl gene ID, so the small number of non-coding entries without a clean
protein mapping (e.g., TERC) carry conservation and interaction annotation
unevenly. Ortholog calls, expression values, and interaction scores inherit the
assumptions and versions of their upstream resources (Ensembl BioMart, GTEx v8,
HPA, STRING v12); because the pipeline pulls from live APIs, re-running it against
newer releases is a single command and will naturally refresh these layers.

The entire database is reproducible end to end from public APIs via `main.sh`,
with figure code preserved in an executable notebook. Natural extensions include
additional model organisms, cell-type-resolved expression from single-cell
atlases, and directional/regulatory edges to complement the undirected STRING
network. As released, the resource provides an integrated, provenance-tracked,
functionally annotated starting point for computational senescence research.

## Figures

- **Figure 1.** Source overlap across the four seed lists (UpSet-style
  intersection sizes with source-membership matrix). — `figures/source_overlap.png`
- **Figure 2.** GTEx v8 median expression heatmap for genes supported by ≥3
  sources, hierarchically clustered on both axes. — `figures/expression_heatmap.png`
- **Figure 3.** STRING v12 high-confidence interaction network (top 60 hubs; node
  size ∝ degree, color = number of confirming sources). — `figures/interaction_network.png`
- **Figure 4.** Net-new annotation coverage by layer and redundancy removed by
  unification. — `figures/coverage_gain.png`
- **Figure 5.** Multi-evidence priority score for the top 20 senescence gene
  candidates. Stacked horizontal bars show the three score components: source
  confirmation (red), network centrality (blue), and evolutionary conservation
  (green). Gene labels are coloured by the number of confirming sources. Scores
  are computed from the unified database without any external data. —
  `figures/priority_score.png`
- **Figure 6.** GTEx v8 median expression (log2 TPM + 1) of the top 15
  priority-scored senescence genes across 15 aging-relevant tissues, with
  independent hierarchical clustering of genes (rows) and tissues (columns). The
  left colour strip encodes the number of confirming sources per gene. —
  `figures/aging_tissue_landscape.png`

## Data and code availability

All code, the reproduction script (`main.sh`), the figure notebook
(`notebooks/figures.ipynb`), and the output database (`results/`) are included in
this repository. Every dataset is fetched from its public API at run time.

## References (key resources)

- Avelar et al. CellAge / Human Ageing Genomic Resources (HAGR). https://genomics.senescence.info/cells/
- Tacutu et al. GenAge / HAGR. https://genomics.senescence.info/genes/
- Saul et al. (2022) SenMayo (SAUL_SEN_MAYO, MSigDB M45803). *Nature Communications*.
- Reactome Cellular Senescence pathway R-HSA-2559583. https://reactome.org
- MyGene.info; Ensembl BioMart; GTEx v8; Human Protein Atlas; STRING v12.
- Acosta JC et al. (2013) CXCR2 and its ligands mediate the senescence-reinforcing autocrine loop. *Nature Cell Biology*, 15, 978–990.
- Birch J & Gil J. (2020) Senescence and the SASP: many therapeutic avenues. *Genes & Development*, 34, 1565–1576.
- Campisi J & d'Adda di Fagagna F. (2007) Cellular senescence: when bad things happen to good cells. *Nature Reviews Molecular Cell Biology*, 8, 729–740.
- Childs BG et al. (2017) Cellular senescence in aging and age-related disease: from mechanisms to therapy. *Nature Medicine*, 23, 1424–1435.
- Clevers H & Nusse R. (2012) Wnt/β-Catenin signaling and disease. *Cell*, 149, 1192–1205.
- Coppe JP et al. (2008) Senescence-associated secretory phenotypes reveal cell-nonautonomous functions of oncogenic RAS and the p53 tumor suppressor. *PLOS Biology*, 6, e301.
- Herranz N & Gil J. (2018) Mechanisms and functions of cellular senescence. *Journal of Clinical Investigation*, 128, 1238–1246.
- Kirkland JL & Tchkonia T. (2017) Cellular senescence: a translational perspective. *EBioMedicine*, 21, 21–28.
- López-Otín C et al. (2013) The hallmarks of aging. *Cell*, 153, 1194–1217.
- Lopes-Paciencia S et al. (2019) The senescence-associated secretory phenotype and its regulation. *Cytokine*, 117, 15–22.
- Rufini A et al. (2013) Senescence and aging: the critical roles of p53. *Oncogene*, 32, 5129–5143.
- Vousden KH & Lane DP. (2007) p53 in health and disease. *Nature Reviews Molecular Cell Biology*, 8, 275–283.
- Yu H et al. (2009) STATs in cancer inflammation and immunity: a leading role for STAT3. *Nature Reviews Cancer*, 9, 798–809.
