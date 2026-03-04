# Cyclic Peptide Binder Design for a Protein Dimer Interface

Computational pipeline for designing cyclic peptide binders targeting the
protein dimer interface of PDB structure 8ZCS. The pipeline spans four
sequential phases: structure preparation, binding pocket prediction,
AlphaFold2-based peptide design, and molecular dynamics validation.

## Pipeline Overview

```
01_structure_prep/  →  02_pocket_analysis/  →  03_design/  →  04_simulation/
   (PDB cleanup)       (AF2BIND hotspots)     (ColabDesign)    (AMBER MD + MM/GBSA)
```

### Phase 1 — Structure Preparation (`01_structure_prep/`)
Input PDB (8zcs.pdb) of the target protein complex.

### Phase 2 — Pocket Analysis (`02_pocket_analysis/`)
Runs AF2BIND to identify binding hotspot residues on the dimer interface.

### Phase 3 — Peptide Design (`03_design/`)
Uses ColabDesign's AlphaFold2 multimer model with custom loss functions to
generate 14-residue cyclic peptide binders.

Post-design analysis pipeline (`03_design/analysis/`):
0. `00_sort_results.py` — parse run logs into `results.csv`, ranked by pLDDT
1. `01_filter_best.py` — filter by hard constraint & pLDDT threshold, copy PDBs
2. `02_pdb_clean.py` — add TER lines, assign chains, run pdb4amber
3. `03_filter_cys.py` — remove designs with single cysteines or unknown residues
4. `04_bsa.py` — buried surface area calculation
5. `05_mpnn.py` — MPNN sequence validation, PSSM correlation, identity metrics
6. `06_merge.py` — merge MPNN results into `results.csv`
7. `07_correlation.py` — statistical correlations

Each step reads and updates a single `results.csv` file, adding columns as the pipeline progresses.

### Phase 4 — Simulation (`04_simulation/`)
AMBER molecular dynamics with MM/GBSA binding energy analysis.

## Installation

### 1. Create the conda environment

```bash
conda env create -f environment.yml
conda activate peptide_design
```

### 2. Install ColabDesign

ColabDesign requires a manual installation. Follow the instructions at
https://github.com/sokrypton/ColabDesign.

### 3. AMBER / AmberTools

Molecular dynamics simulations in Phase 4 require AMBER (pmemd.cuda) and
AmberTools (tleap, pdb4amber, MMPBSA.py). See https://ambermd.org/.

### 4. AlphaFold2 weights

Download AlphaFold2 model parameters and update the path in `config.yaml`.

## Configuration

All user-configurable paths and design parameters are in **`config.yaml`** at the
repository root. Edit this file before running any pipeline step:

```yaml
paths:
  af2_params_dir: "/path/to/alphafold/weights"
  af2bind_params: "/path/to/af2bind_params/attempt_7_2k_lam0-03"
  target_pdb: "01_structure_prep/8zcs.pdb"
```

## Usage

### Design phase

```bash
# Submit via SLURM
sbatch 03_design/submit.sh

# Or run directly
python3 -u 03_design/run_design.py -n <num_runs> -o <output_dir>
```

### Analysis pipeline (run sequentially)

Step 05 requires a GPU and the full `peptide_design` environment with ColabDesign. All other steps require only the base scientific Python stack.

```bash
cd 03_design/analysis
python3 00_sort_results.py        # parse design run logs → results.csv
python3 01_filter_best.py         # filter by pLDDT + hard constraint
python3 02_pdb_clean.py           # clean PDBs via pdb4amber
python3 03_filter_cys.py          # remove designs with Cys/unknown residues
python3 04_bsa.py                 # buried surface area (MDTraj)
python3 05_mpnn.py                # MPNN sequence validation (ColabDesign, GPU)
python3 06_merge.py               # merge MPNN results into results.csv
# python3 07_correlation.py       # optional: statistical correlations (scipy)
```

### Simulation setup and execution

```bash
cd 04_simulation/<batch>/
bash ../setup_scripts/create_dir.sh -a 1-153 -b 154-303 -p 304-317
python3 ../setup_scripts/new_mask.py <range>
bash ../setup_scripts/leap.sh
bash ../setup_scripts/submit_all.sh
```

## License

MIT
