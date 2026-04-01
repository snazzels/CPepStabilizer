#!/bin/bash
#SBATCH --job-name=mmgbsa-BL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:0
#SBATCH --time=1:00:00
#SBATCH --mem=10000MB
module load ambertools

echo "===INFORMATION==="
echo "Script dir:  $(dirname $0)"
echo "Current dir: $PWD"
echo "Argument:    $1"
echo "Amberhome:   $AMBERHOME"
echo "Python:      `type python`"
echo "CUDA home:   $CUDA_HOME"
echo "Nodes:       $SLURM_JOB_NODELIST"
echo "======"

mpirun -np 8 MMPBSA.py.MPI -i ../gb.in -y ../LB.nc -cp ../LB.prmtop -sp ../LB.prmtop -rp ../B.prmtop -lp ../L.prmtop -o BL_output.dat -do BL_decomp.csv

