#!/usr/bin/env bash
set -euo pipefail

# Read interpreter paths from config.yaml (can be overridden via env vars)
_cfg() { python3 -c "
import re, sys
key, path = sys.argv[1], sys.argv[2]
in_env = False
for line in open(path):
    if re.match(r'^environments:', line): in_env = True
    elif in_env and re.match(r'^\S', line): in_env = False
    elif in_env and re.match(rf'\s+{key}:', line):
        print(line.split(':', 1)[1].strip().strip('\"'))
        sys.exit()
" "$1" "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/config.yaml"; }

PYTHON="${PYTHON:-$(_cfg python)}"
PYTHON_GPU="${PYTHON_GPU:-$(_cfg python_gpu)}"

# Prepend GPU env's bin to PATH so JAX finds the bundled ptxas (conda ships 12.x,
# system ptxas may be too old for newer jaxlib builds)
GPU_ENV_BIN="$(dirname "$PYTHON_GPU")"
export PATH="$GPU_ENV_BIN:$PATH"
RUNS_DIR="../../example_design_runs"
SELECT_ONLY=0

while [[ $# -gt 0 ]]; do
  case $1 in
    --runs_dir) RUNS_DIR="$2"; shift 2 ;;
    --select) SELECT_ONLY=1; shift ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

cd "$(dirname "${BASH_SOURCE[0]}")"

if [[ $SELECT_ONLY -eq 1 ]]; then
  echo "[06] Selecting designs for simulation..."
  "$PYTHON" 06_select_for_sim.py
  echo "=== Selection complete. ==="
  exit 0
fi

echo "[00] Parsing run logs and filtering by pLDDT..."
"$PYTHON" 00_sort_results.py --runs_dir "$RUNS_DIR"
echo "[01] Cleaning PDBs via pdb4amber..."
"$PYTHON" 01_pdb_clean.py
echo "[02] Filtering Cys/unknown residues..."
"$PYTHON" 02_filter_cys.py
echo "[03] Calculating buried surface area..."
"$PYTHON" 03_bsa.py
echo "[04] Running MPNN sequence validation (GPU)..."
"$PYTHON_GPU" 04_mpnn.py
echo "[Summary] Generating plots and table..."
"$PYTHON" 05_summary.py
echo "=== Pipeline complete. ==="
echo "Inspect summary_plots.png, set thresholds in config.yaml, then run:"
echo "  bash run_pipeline.sh --select"
