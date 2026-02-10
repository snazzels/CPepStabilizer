#!/bin/bash
#SBATCH --job-name=test_name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --exclude=t38cn030,t38cn038,t38cn028,t38cn015,t38cn024,t38cn027,t38cn016,t38cn018,t38cn031,t38cn020,t38cn021,t38cn032,t38cn019,t38cn022,t38cn054,t38cn053
#SBATCH --time=1-23:00:00

module load cuda
module load pmemd
module load ambertools

export CUDA_VISIBLE_DEVICES=0
# export CUDA_LAUNCH_BLOCKING=1


echo "===INFORMATION==="
echo "Script dir:  $(dirname $0)"
echo "Current dir: $PWD"
echo "Argument:    $1"
echo "Amberhome:   $AMBERHOME"
echo "Python:      `type python`"
echo "CUDA home:   $CUDA_HOME"
echo "GPU Info:    $(nvidia-smi -L)"
echo "Nodes:       $SLURM_JOB_NODELIST"
echo "======"


echo "==========Simulation START============="
echo "Date of simulation start: $(date)"



# Run Simulation
bash ./run.sh


# MMGBSA
cd ../mmgbsa
bash ./single.sh -a 1-210 -b 277-376 -p 377-384
#bash ./submit.sh

#bash ./analyze_mmgbsa.py

echo "==========Simulation DONE=============="
echo "Date of completion: $(date)"
