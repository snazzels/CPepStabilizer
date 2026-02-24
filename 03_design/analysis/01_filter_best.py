#!/usr/bin/env python3
"""
Step 01: Filter designs by hard constraint and pLDDT threshold, copy PDB files.

Reads results.csv, filters rows, adds pdb_filename column, overwrites results.csv.
"""

import os
import shutil
import pandas as pd
import yaml
from pathlib import Path

_repo_root = Path(__file__).resolve().parents[2]
with open(_repo_root / "config.yaml") as _f:
    _config = yaml.safe_load(_f)


def process_design_results(input_csv, best_top_dir, design_runs_dir):
    plddt_threshold = _config["analysis"]["plddt_threshold"]

    os.makedirs(best_top_dir, exist_ok=True)

    df = pd.read_csv(input_csv)
    filtered_df = df[(df['hard'] == 1) & (df['plddt'] >= plddt_threshold)].copy()

    filenames = []
    for _, row in filtered_df.iterrows():
        run_id = row['run_id']
        model_number = row['models'] - 1
        trajectory_dir = os.path.join(design_runs_dir, run_id, 'trajectory')
        pdb_file = os.path.join(trajectory_dir, f'model_step_{model_number}.pdb')

        if os.path.isfile(pdb_file):
            new_filename = f'{run_id}_model_step_{model_number}.pdb'
            shutil.copy(pdb_file, os.path.join(best_top_dir, new_filename))
            filenames.append(new_filename)
            print(f'Copied {pdb_file} to {best_top_dir}')
        else:
            filenames.append(None)
            print(f'File not found: {pdb_file}')

    filtered_df['pdb_filename'] = filenames
    filtered_df.to_csv(input_csv, index=False)

    print(f'\n{len(filtered_df)} designs passed filter (pLDDT >= {plddt_threshold}, hard == 1)')
    print(f'Updated {input_csv}')
    return filtered_df


def main():
    result_df = process_design_results(
        input_csv='results.csv',
        best_top_dir='best_top',
        design_runs_dir='../design_runs',
    )
    print("\nFiltered Results:")
    print(result_df)


if __name__ == "__main__":
    main()
