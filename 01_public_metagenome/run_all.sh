#!/usr/bin/env bash
set -euo pipefail

echo "[1/2] Main feature selection"
python scripts/01_main_stage_feature_selection.py

echo "[2/2] G. morbillorum stage-wise detection rate and relative abundance"
python scripts/02_gm_stage_detection_abundance.py

echo "All manuscript main analyses completed."
