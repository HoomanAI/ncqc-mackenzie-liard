#!/usr/bin/env bash
# Runs the full NCQC workflow from the raw workbooks. Takes roughly 10-15 minutes.
set -euo pipefail
cd "$(dirname "$0")"
python src/01_load_data.py
python src/03_error_injection.py
python src/04_make_figures.py
python src/05_sensitivity_and_baselines.py
echo "Done. Outputs are in data/processed/ and figures/."
