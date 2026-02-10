# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Computational biochemistry research pipeline for designing cyclic peptide binders targeting a protein dimer interface (PDB: 8ZCS). The pipeline spans four sequential phases: structure preparation, binding pocket prediction, AlphaFold2-based peptide design, and molecular dynamics validation.

## Configuration

All user-configurable paths and design parameters live in **`config.yaml`** at the repo root. Scripts load this file at startup to resolve paths relative to the repository root.

## Pipeline Architecture

The project follows a strict four-phase pipeline where each phase's outputs feed into the next:

```
01_structure_prep/  →  02_pocket_analysis/  →  03_design/  →  04_simulation/
   (PDB cleanup)       (AF2BIND hotspots)     (ColabDesign)    (AMBER MD + MM/GBSA)
```

### Phase 1: Structure Preparation (`01_structure_prep/`)
Input PDB (8zcs.pdb) of the target protein complex.

### Phase 2: Pocket Analysis (`02_pocket_analysis/`)
`af2bind.py` runs AF2BIND inference to identify binding hotspot residues. Paths loaded from `config.yaml`. Outputs `results.csv` with per-residue binding probabilities.

### Phase 3: Peptide Design (`03_design/`)
Core design engine. `run_design.py` orchestrates ColabDesign's AlphaFold2 multimer model with custom loss functions to generate 14-residue cyclic peptide binders. All design parameters are loaded from `config.yaml`.

**Custom loss functions** (`functions_design/`):
- `constants.py`: Shared Van der Waals radii dictionary
- `loss_functions.py`: `PeptideLoss` class — cyclic offset, COM distance from hotspot, cis/trans penalty, dihedral computation
- `sasa_functions.py`: Vectorized Shrake-Rupley SASA via JAX vmap, 3-component BSA loss (weights: w1=1.0, w2=0.5, w3=2.0)

**Post-design analysis** (`03_design/analysis/`) runs as a numbered sequential pipeline:
1. `01_filter_best.py` — filter by hard constraint & pLDDT threshold
2. `02_pdb_clean.py` — add TER lines, assign chains, run pdb4amber
3. `03_filter_cys.py` — remove designs with single cysteines or unknown residues
4. `04_bsa.py` — buried surface area calculation
5. `05_mpnn.py` — MPNN sequence validation, PSSM correlation, identity metrics
6. `06_merge.py` — merge design metrics with MPNN results
7. `07_correlation.py` — statistical correlations

Additional analysis scripts: `composition.py`, `logo.py`, `perplexity.py`, `paper_filter.py`, `paper_plot.py`.

### Phase 4: Simulation (`04_simulation/`)
AMBER molecular dynamics with MM/GBSA binding energy analysis. 3 replicate runs per batch, multiple design models per run.

**Setup pipeline** (`setup_scripts/`):
- `create_dir.sh -a <chainA_range> -b <chainB_range> -p <peptide_range>` — template expansion for simulation directories
- `new_mask.py` — residue mask generation for restraints
- `leap.sh` — LEaP topology builder (solvation with 0.15M NaCl)
- `submit_all.sh` — batch SLURM submission

**Analysis** (per batch directory):
- `average.py` — MM/GBSA averaging across 3 replicas
- `coop.py` — cooperativity analysis

## Running Commands

### Design phase
```bash
# Submit via SLURM
sbatch 03_design/submit.sh

# Direct run
python3 -u 03_design/run_design.py -n <num_runs> -o <output_dir>
```

### Analysis pipeline (sequential, run in order)
```bash
cd 03_design/analysis
python3 01_filter_best.py
python3 02_pdb_clean.py
python3 03_filter_cys.py
python3 04_bsa.py
python3 05_mpnn.py
python3 06_merge.py
python3 07_correlation.py
```

### Simulation setup and execution
```bash
cd 04_simulation/<batch>/
bash ../setup_scripts/create_dir.sh -a 1-153 -b 154-303 -p 304-317
python3 ../setup_scripts/new_mask.py <range>
bash ../setup_scripts/leap.sh
bash ../setup_scripts/submit_all.sh
```

## Environment

A single conda environment covers all dependencies:

```bash
conda env create -f environment.yml
conda activate peptide_design
```

Additionally requires:
- **ColabDesign**: AF2 design framework (custom install)
- **AMBER suite**: pdb4amber, tleap, pmemd.cuda, MMPBSA.py
- **SLURM**: Job scheduler for HPC cluster

## Key Dependencies

- **ColabDesign**: AF2 design framework (custom install, not pip)
- **JAX/jaxlib**: JIT compilation and autodiff for loss functions
- **AMBER suite**: pdb4amber, tleap, mmpbsa_py for MD
- **BioPython, MDTraj, Pandas, NumPy, SciPy, Matplotlib, PyYAML**

## Important Notes

- No automated tests or CI/CD — validation is manual via CSV inspection and structure visualization
- Design runs are stored in timestamped directories: `run_NNN_YYYYMMDD_HHMMSS/`
- The `LogRedirector` class in `run_design.py` captures stdout to per-run `run.log` files
- Loss functions use JAX tensors and are differentiable for AF2 backpropagation
