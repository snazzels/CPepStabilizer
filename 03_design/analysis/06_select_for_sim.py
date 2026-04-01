#!/usr/bin/env python3
"""
Step 09: Select designs for simulation based on seq-identity and min-BSA thresholds.

Thresholds are read from config.yaml (analysis.seq_id_threshold and
analysis.bsa_min_threshold). pLDDT >= 0.8 is already enforced by step 01.

Passing PDB files are copied from pdbs/ to sim_pdbs/.
"""

import os
import shutil
import pandas as pd
from _load_config import load_config

_config = load_config()
seq_id_threshold = _config["analysis"]["seq_id_threshold"]
bsa_min_threshold = _config["analysis"]["bsa_min_threshold"]

df = pd.read_csv("results.csv")

if "bsa_a" in df.columns and "bsa_b" in df.columns:
    df["bsa_min"] = df[["bsa_a", "bsa_b"]].min(axis=1)

missing = [c for c in ("Avg Exact Identity", "bsa_min") if c not in df.columns]
if missing:
    raise SystemExit(f"results.csv is missing columns: {missing}. Run steps 04–06 first.")

mask = (
    (df["Avg Exact Identity"] >= seq_id_threshold) &
    (df["bsa_min"] >= bsa_min_threshold)
)
selected = df[mask].copy()

os.makedirs("sim_pdbs", exist_ok=True)
for fname in selected["pdb_filename"].dropna():
    src = os.path.join("pdbs", fname)
    if os.path.isfile(src):
        shutil.copy(src, os.path.join("sim_pdbs", fname))
    else:
        print(f"  Warning: {src} not found, skipping.")

print(f"Thresholds: seq_id >= {seq_id_threshold},  min(BSA) >= {bsa_min_threshold} Å²")
print(f"Selected {len(selected)} / {len(df)} designs → sim_pdbs/")
print()
print("Next steps — from the repo root:")
print("  cd 04_simulation/run_sim")
print("  bash setup.sh")
print("  bash ../setup_scripts/submit_all.sh")
