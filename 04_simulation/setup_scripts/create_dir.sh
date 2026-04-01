#!/bin/bash
set -euo pipefail

pdb_dir="../../03_design/analysis/sim_pdbs/"
sim_dir="../sim_scripts"

# Parse optional --pdb_dir override
while [[ $# -gt 0 ]]; do
  case $1 in
    --pdb_dir) pdb_dir="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

# Extract residue ranges for chains A, B, C from a PDB file.
# Prints three lines:  A start end  /  B start end  /  C start end
_ranges() {
  python3 - "$1" <<'EOF'
import sys
chains = {}
for line in open(sys.argv[1]):
    if line.startswith("ATOM"):
        ch  = line[21]
        res = int(line[22:26])
        if ch not in chains:
            chains[ch] = [res, res]
        else:
            if res < chains[ch][0]: chains[ch][0] = res
            if res > chains[ch][1]: chains[ch][1] = res
for ch in ("A", "B", "C"):
    lo, hi = chains[ch]
    print(f"{ch} {lo} {hi}")
EOF
}

for pdb_file in "$pdb_dir"/*.pdb; do
    pdb_filename=$(basename "$pdb_file" .pdb)
    sim_subdir="$pdb_filename"

    # Derive peptide residue range from this PDB (needed for cyclic bond in tleap)
    ranges=$(_ranges "$pdb_file")
    p_range=$(echo "$ranges" | awk '$1=="C"{print $2"-"$3}')
    p1=$(echo "$p_range" | cut -d- -f1)
    p2=$(echo "$p_range" | cut -d- -f2)

    echo "$pdb_filename  pep:$p_range"

    cp -r "$sim_dir" "$sim_subdir"
    cp "$pdb_file" "$sim_subdir/top_prep/"

    # tleap: load correct PDB and set cyclic bond
    sed -i "s/mol=loadpdb run_1.pdb/mol=loadpdb $pdb_filename.pdb/" \
        "$sim_subdir/top_prep/tleap.in"
    sed -i "s/bond mol\.[0-9]*\.N mol\.[0-9]*\.C/bond mol.${p1}.N mol.${p2}.C/" \
        "$sim_subdir/top_prep/tleap.in"

    # SLURM job name
    sed -i "s/#SBATCH --job-name=test_name/#SBATCH --job-name=$pdb_filename/" \
        "$sim_subdir/mdin/submit.sh"

done
