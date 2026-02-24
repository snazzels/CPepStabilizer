#!/usr/bin/env python3
"""
Step 06: Merge MPNN sequence identity results into results.csv.

Joins on pdb_filename and overwrites results.csv with the combined metrics.
"""

import pandas as pd

df = pd.read_csv("results.csv")
mpnn_df = pd.read_csv("output/sequence_identity_results.csv")

merged_df = pd.merge(df, mpnn_df, on="pdb_filename", how="inner")
merged_df.to_csv("results.csv", index=False)

print(f"Merged {len(merged_df)} rows into results.csv")
