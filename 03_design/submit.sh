#!/bin/bash
#SBATCH --job-name=peptide_design

cd "$(dirname "$0")"

python3 -u run_design.py -n 160
