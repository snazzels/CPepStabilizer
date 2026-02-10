#!/usr/bin/env python3

import pandas as pd

# Paths and directories
file1_path = "filtered_file_bsa.csv"
file2_path = "output/sequence_identity_results.csv"

# Load the CSV files
df1 = pd.read_csv(file1_path)
df2 = pd.read_csv(file2_path)

# Standardize filenames by removing .pdb and .csv extensions
df1["pdb_filename"] = df1["pdb_filename"].str.replace(".pdb", "", regex=False)
df2["pdb_filename"] = df2["File"].str.replace(".csv", "", regex=False)

# Merge based on the cleaned pdb_filename column
merged_df = pd.merge(df1, df2, on="pdb_filename", how="inner")

# Drop the 'File' column
merged_df = merged_df.drop(columns=['File'])

# Save the merged results
merged_df.to_csv("sasa.csv", index=False)
