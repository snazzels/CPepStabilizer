#!/usr/bin/env python3
"""
Step 03: Filter out designs containing single cysteines (C) or unknown residues (X).

Reads PDB files from pdbs/, deletes offending files, logs them to cys_filtered.csv,
and removes the corresponding rows from results.csv.
"""

import os
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


def filter_cys(pdbs_dir="pdbs", csv_path="results.csv", log_path="cys_filtered.csv"):
    """Filter PDBs with C or X in chain C, delete them, log to CSV, and update results.csv."""
    df = pd.read_csv(csv_path)
    removed_rows = []

    pdb_files = [f for f in os.listdir(pdbs_dir) if f.lower().endswith(".pdb")]

    for pdb_file in pdb_files:
        pdb_path = os.path.join(pdbs_dir, pdb_file)
        sequence = get_chain_c_sequence(pdb_path)

        if "C" in sequence or "X" in sequence:
            os.remove(pdb_path)
            removed_rows.append({"pdb_filename": pdb_file, "sequence": sequence, "reason": "C or X in peptide"})
            print(f"Deleted {pdb_file} (sequence: {sequence})")
        else:
            print(f"OK: {pdb_file} (sequence: {sequence})")

    if removed_rows:
        removed_files = {r["pdb_filename"] for r in removed_rows}
        pd.DataFrame(removed_rows).to_csv(log_path, index=False)
        initial_len = len(df)
        df = df[~df["pdb_filename"].isin(removed_files)]
        df.to_csv(csv_path, index=False)
        print(f"\nRemoved {initial_len - len(df)} rows from {csv_path}, logged to {log_path}")
    else:
        print("\nNo designs with C or X found.")

    print(f"Remaining designs: {len(df)}")


if __name__ == "__main__":
    filter_cys(pdbs_dir="pdbs")
