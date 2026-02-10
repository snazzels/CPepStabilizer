#!/usr/bin/env python3
"""
Step 03: Filter out designs containing single cysteines (C) or unknown residues (X).

Reads PDB files from cleaned_top/, moves offending files to CYS/, and updates
filtered_best.csv by removing the corresponding rows.
"""

import os
import shutil
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


def get_chain_c_sequence(pdb_file):
    """Extract the amino acid sequence of chain C (peptide) from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("PDB", pdb_file)
    sequence = ""
    first_model = next(structure.get_models())
    for chain in first_model:
        if chain.id == "C":
            for residue in chain:
                if residue.id[0] == " ":  # Exclude heteroatoms
                    sequence += seq1(residue.get_resname())
    return sequence


def filter_cys(cleaned_top_dir="cleaned_top", csv_path="filtered_best.csv"):
    """Filter PDBs with C or X in chain C and update the CSV accordingly."""
    cys_dir = os.path.join(os.path.dirname(cleaned_top_dir.rstrip("/")), "CYS")
    os.makedirs(cys_dir, exist_ok=True)

    df = pd.read_csv(csv_path)
    removed_files = set()

    pdb_files = [f for f in os.listdir(cleaned_top_dir) if f.lower().endswith(".pdb")]

    for pdb_file in pdb_files:
        pdb_path = os.path.join(cleaned_top_dir, pdb_file)
        sequence = get_chain_c_sequence(pdb_path)

        if "C" in sequence or "X" in sequence:
            dest_path = os.path.join(cys_dir, pdb_file)
            shutil.move(pdb_path, dest_path)
            removed_files.add(pdb_file)
            print(f"Moved {pdb_file} to {cys_dir} (sequence: {sequence})")
        else:
            print(f"OK: {pdb_file} (sequence: {sequence})")

    # Update CSV: remove rows whose PDB was moved
    if removed_files:
        initial_len = len(df)
        df = df[~df["pdb_filename"].isin(removed_files)]
        df.to_csv(csv_path, index=False)
        print(f"\nRemoved {initial_len - len(df)} rows from {csv_path}")
    else:
        print("\nNo designs with C or X found.")

    print(f"Remaining designs: {len(df)}")


if __name__ == "__main__":
    filter_cys()
