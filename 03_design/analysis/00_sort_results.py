#!/usr/bin/env python3
"""
Step 00: Parse design run logs, filter by pLDDT threshold, and copy PDB files.

Reads run.log files from runs_dir, filters by hard == 1 and pLDDT >= threshold
(from config.yaml), copies passing PDB files to pdbs/, and writes results.csv.
"""

import os
import re
import shutil
import argparse
import pandas as pd
from _load_config import load_config

_config = load_config()


def extract_data_from_log(log_path):
    data = []
    run_id = os.path.basename(os.path.dirname(log_path))
    with open(log_path, 'r') as f:
        for line in f:
            match = re.match(
                r"(\d+) models \[(\d+)\] recycles (\d+) hard (\d+) soft ([\d.]+) "
                r"temp ([\d.]+) loss ([\d.]+) i_con_1 ([\d.]+) i_con_2 ([\d.]+) "
                r"plddt ([\d.]+) ptm ([\d.]+) i_ptm ([\d.]+)", line)
            if match:
                data.append({
                    "run_id": run_id,
                    "models": int(match.group(1)),
                    "index": int(match.group(2)),
                    "recycles": int(match.group(3)),
                    "hard": int(match.group(4)),
                    "soft": float(match.group(5)),
                    "temp": float(match.group(6)),
                    "loss": float(match.group(7)),
                    "i_con_1": float(match.group(8)),
                    "i_con_2": float(match.group(9)),
                    "plddt": float(match.group(10)),
                    "ptm": float(match.group(11)),
                    "i_ptm": float(match.group(12)),
                })
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--runs_dir', default='../design_runs')
    parser.add_argument('--output', default='results.csv')
    args = parser.parse_args()

    plddt_threshold = _config["analysis"]["plddt_threshold"]
    pdbs_dir = "pdbs"
    os.makedirs(pdbs_dir, exist_ok=True)

    all_data = []
    for run_dir in sorted(os.listdir(args.runs_dir)):
        log_path = os.path.join(args.runs_dir, run_dir, "run.log")
        if os.path.isfile(log_path):
            all_data.extend(extract_data_from_log(log_path))

    all_data.sort(key=lambda x: x['plddt'], reverse=True)

    df = pd.DataFrame(all_data)
    filtered = df[(df['hard'] == 1) & (df['plddt'] >= plddt_threshold)].copy()

    filenames = []
    for _, row in filtered.iterrows():
        model_number = row['models'] - 1
        pdb_src = os.path.join(args.runs_dir, row['run_id'], 'trajectory',
                               f'model_step_{model_number}.pdb')
        if os.path.isfile(pdb_src):
            new_filename = f"{row['run_id']}_model_step_{model_number}.pdb"
            shutil.copy(pdb_src, os.path.join(pdbs_dir, new_filename))
            filenames.append(new_filename)
        else:
            filenames.append(None)
            print(f'File not found: {pdb_src}')

    filtered['pdb_filename'] = filenames
    filtered.to_csv(args.output, index=False)

    print(f"Parsed {len(all_data)} entries from logs")
    print(f"{len(filtered)} designs passed filter (pLDDT >= {plddt_threshold}, hard == 1) → {pdbs_dir}/")


if __name__ == "__main__":
    main()
