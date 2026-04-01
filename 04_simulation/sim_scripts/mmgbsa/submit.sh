#!/bin/bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

sbatch --chdir="$DIR/AL_out"   "$DIR/AL_out/AL.sh"
sbatch --chdir="$DIR/BL_out"   "$DIR/BL_out/BL.sh"
sbatch --chdir="$DIR/AB_L_out" "$DIR/AB_L_out/AB_L.sh"
sbatch --chdir="$DIR/AB_out"   "$DIR/AB_out/AB.sh"
sbatch --chdir="$DIR/AL_B_out" "$DIR/AL_B_out/AL_B.sh"
