#!/usr/bin/env bash
# One-click reproduction of the unified senescence/SASP gene database.
#   ./main.sh            run the full data pipeline + regenerate all figures
#   ./main.sh --no-figs  data pipeline only
set -euo pipefail
cd "$(dirname "$0")"

PY="${PYTHON:-python3}"

echo ">> Installing Python dependencies"
$PY -m pip install -q -r requirements.txt

echo ">> Running data pipeline (steps 1-8)"
$PY run_pipeline.py

if [[ "${1:-}" != "--no-figs" ]]; then
  echo ">> Regenerating figures from notebooks/figures.ipynb"
  $PY -m jupyter nbconvert --to notebook --execute --inplace \
      --ExecutePreprocessor.timeout=600 notebooks/figures.ipynb
  echo ">> Figures written to figures/"
fi

echo ">> DONE. Outputs:"
echo "   results/senescence_sasp_database.json   (main database)"
echo "   results/senescence_sasp_database.csv"
echo "   results/summary_report.md"
echo "   results/build_metrics.json"
echo "   figures/*.png"
