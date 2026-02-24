#!/usr/bin/env python3
"""
Step 04: Calculate buried surface area (BSA) for each design.

Reads results.csv, adds bsa_a, bsa_b, sasa_peptide columns, overwrites results.csv.
"""

import os
import mdtraj as md
import numpy as np
import pandas as pd
from _load_config import load_config

_config = load_config()

probe_radius = _config["sasa"]["probe_radius"]
n_sphere_points = _config["sasa"]["n_sphere_points"]

df = pd.read_csv("results.csv")

bsa_a_values = []
bsa_b_values = []
sasa_peptide_in_complex_values = []

best_top_dir = 'cleaned_top'

for _, row in df.iterrows():
    pdb_filename = row['pdb_filename']

    if pd.isna(pdb_filename):
        bsa_a_values.append(None)
        bsa_b_values.append(None)
        sasa_peptide_in_complex_values.append(None)
        continue

    pdb_file_path = os.path.join(best_top_dir, pdb_filename)

    if not os.path.isfile(pdb_file_path):
        print(f"File {pdb_file_path} not found, skipping.")
        bsa_a_values.append(None)
        bsa_b_values.append(None)
        sasa_peptide_in_complex_values.append(None)
        continue

    print(f"Processing file: {pdb_filename}")

    traj = md.load(pdb_file_path)

    chain_a_atoms = traj.top.select("chainid 0")
    chain_b_atoms = traj.top.select("chainid 1")
    peptide_atoms = traj.top.select("chainid 2")
    chain_a_peptide_atoms = traj.top.select("chainid 0 or chainid 2")
    chain_b_peptide_atoms = traj.top.select("chainid 1 or chainid 2")

    chain_a_traj = traj.atom_slice(chain_a_atoms)
    chain_b_traj = traj.atom_slice(chain_b_atoms)
    peptide_traj = traj.atom_slice(peptide_atoms)
    chain_a_peptide_traj = traj.atom_slice(chain_a_peptide_atoms)
    chain_b_peptide_traj = traj.atom_slice(chain_b_peptide_atoms)

    sasa_chain_a = md.shrake_rupley(chain_a_traj, probe_radius=probe_radius, n_sphere_points=n_sphere_points).sum(axis=1) * 100
    sasa_chain_b = md.shrake_rupley(chain_b_traj, probe_radius=probe_radius, n_sphere_points=n_sphere_points).sum(axis=1) * 100
    sasa_peptide_isolated = md.shrake_rupley(peptide_traj, probe_radius=probe_radius, n_sphere_points=n_sphere_points).sum(axis=1) * 100
    sasa_chain_a_peptide = md.shrake_rupley(chain_a_peptide_traj, probe_radius=probe_radius, n_sphere_points=n_sphere_points).sum(axis=1) * 100
    sasa_chain_b_peptide = md.shrake_rupley(chain_b_peptide_traj, probe_radius=probe_radius, n_sphere_points=n_sphere_points).sum(axis=1) * 100

    bsa_a = sasa_chain_a + sasa_peptide_isolated - sasa_chain_a_peptide
    bsa_b = sasa_chain_b + sasa_peptide_isolated - sasa_chain_b_peptide

    sasa_total = md.shrake_rupley(traj, probe_radius=probe_radius, n_sphere_points=n_sphere_points)
    sasa_peptide_in_complex = sasa_total[:, peptide_atoms].sum(axis=1) * 100

    bsa_a_values.append(bsa_a[0])
    bsa_b_values.append(bsa_b[0])
    sasa_peptide_in_complex_values.append(sasa_peptide_in_complex[0])

    print(f"BSA between Chain A and Peptide: {bsa_a[0]:.2f} Å²")
    print(f"BSA between Chain B and Peptide: {bsa_b[0]:.2f} Å²")
    print(f"SASA of Peptide (in complex): {sasa_peptide_in_complex[0]:.2f} Å²")
    print("-" * 50)

df['bsa_a'] = bsa_a_values
df['bsa_b'] = bsa_b_values
df['sasa_peptide'] = sasa_peptide_in_complex_values

df.to_csv("results.csv", index=False)
print("Updated results.csv with BSA columns")
