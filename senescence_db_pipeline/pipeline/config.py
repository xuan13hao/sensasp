"""Central configuration and shared constants for the senescence/SASP database pipeline."""
from pathlib import Path

# --- Paths (all relative to the project root) ---
ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
DATA_PROC = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
for _d in (DATA_RAW, DATA_PROC, FIGURES, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)

# --- Seed sources ---
SOURCES = ["CellAge", "GenAge", "SenMayo", "Reactome"]

# Source download endpoints (public, no authentication required)
URL_CELLAGE = "https://genomics.senescence.info/cells/cellAge.zip"
URL_GENAGE = "https://genomics.senescence.info/genes/human_genes.zip"
URL_SENMAYO = ("https://www.gsea-msigdb.org/gsea/msigdb/human/"
               "download_geneset.jsp?geneSetName=SAUL_SEN_MAYO&fileType=json")
REACTOME_PATHWAY_ID = "R-HSA-2559583"          # Cellular Senescence
URL_REACTOME_PARTS = f"https://reactome.org/ContentService/data/participants/{REACTOME_PATHWAY_ID}"

# --- Annotation service endpoints ---
URL_MYGENE = "https://mygene.info/v3/query"
URL_GTEX = "https://gtexportal.org/api/v2"
URL_HPA = "https://www.proteinatlas.org/api/search_download.php"
URL_STRING = "https://string-db.org/api"
URL_BIOMART = "https://www.ensembl.org/biomart/martservice"

# --- Parameters ---
GTEX_DATASET = "gtex_v8"
STRING_SPECIES = 9606
STRING_MIN_SCORE = 700            # high-confidence combined score (0-1000)
BIOMART_CHUNK = 200               # Ensembl IDs per BioMart request
MYGENE_CHUNK = 900                # symbols per MyGene POST (<=1000)
REQUEST_TIMEOUT = 60              # seconds
REQUEST_PAUSE = 0.34             # polite delay between per-gene HPA calls (~3 req/s)

# Representative tissues for the expression heatmap panel
HEATMAP_TISSUES = [
    "Cells_Cultured_fibroblasts", "Whole_Blood", "Lung", "Liver", "Kidney_Cortex",
    "Skin_Sun_Exposed_Lower_leg", "Artery_Aorta", "Heart_Left_Ventricle",
    "Brain_Cortex", "Muscle_Skeletal", "Adipose_Subcutaneous", "Spleen",
    "Pancreas", "Colon_Transverse",
]
