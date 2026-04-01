#!/bin/bash
#SBATCH --job-name=mmgbsa-AB_L
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:0
#SBATCH --time=6:00:00
#SBATCH --mem=10000MB
#SBATCH --exclude=t38cn023,t38cn019,t38cn020,t38cn018,t38cn021,t38cn038,t38cn016,t38cn015,t38cn036,t38cn037

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

mpirun -np 8 MMPBSA.py.MPI -i ../gb.in -y ../ALB.nc -cp ../ALB.prmtop -sp ../ALB.prmtop -rp ../AB.prmtop -lp ../L.prmtop -o AB_L_output.dat -do AB_L_decomp.csv

