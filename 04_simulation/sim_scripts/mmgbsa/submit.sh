#! /bin/bash

cd AL_out/
sbatch AL.sh
cd ../BL_out/
sbatch BL.sh
cd ../AB_L_out
sbatch AB_L.sh

