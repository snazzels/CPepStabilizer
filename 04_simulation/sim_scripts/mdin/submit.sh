#!/bin/bash
#SBATCH --job-name=test_name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
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



# Derive chain residue ranges from the PDB in top_prep/
pdb_file=$(ls ../top_prep/*.pdb | head -1)
ranges=$(python3 - "$pdb_file" <<'EOF'
import sys
chains = {}
for line in open(sys.argv[1]):
    if line.startswith("ATOM"):
        ch = line[21]; res = int(line[22:26])
        if ch not in chains: chains[ch] = [res, res]
        else:
            if res < chains[ch][0]: chains[ch][0] = res
            if res > chains[ch][1]: chains[ch][1] = res
for ch in ("A", "B", "C"):
    lo, hi = chains[ch]
    print(f"{ch} {lo} {hi}")
EOF
)
a_range=$(echo "$ranges" | awk '$1=="A"{print $2"-"$3}')
b_range=$(echo "$ranges" | awk '$1=="B"{print $2"-"$3}')
p_range=$(echo "$ranges" | awk '$1=="C"{print $2"-"$3}')

# Run Simulation
bash ./run.sh


# MMGBSA: prepare trajectories then submit AL/BL/AB_L jobs
cd ../mmgbsa
bash ./single.sh -a "$a_range" -b "$b_range" -p "$p_range"
bash ./submit.sh

#bash ./analyze_mmgbsa.py

echo "==========Simulation DONE=============="
echo "Date of completion: $(date)"
