#!/bin/bash

# Default values
pdb_dir="best_top/"
sim_dir="../sim_scripts"
a_range=""
b_range=""
p_range=""

# Function to display help message
usage() {
    echo "Usage: $0 -a startA-endA -b startB-endB -p startP-endP"
    echo "  -a   Specify range for -a in the format start-end (e.g., 1-276)"
    echo "  -b   Specify range for -b in the format start-end (e.g., 277-376)"
    echo "  -p   Specify range for -p in the format start-end (e.g., 377-384)"
    echo "  -h   Display this help message"
    exit 1
}

# Parse command-line options
while getopts ":a:b:p:h" opt; do
    case ${opt} in
        a)
            if [[ $OPTARG =~ ^[0-9]+-[0-9]+$ ]]; then
                a_range="$OPTARG"
            else
                echo "Error: Invalid format for -a. Use start-end (e.g., 1-276)."
                exit 1
            fi
            ;;
        b)
            if [[ $OPTARG =~ ^[0-9]+-[0-9]+$ ]]; then
                b_range="$OPTARG"
            else
                echo "Error: Invalid format for -b. Use start-end (e.g., 277-376)."
                exit 1
            fi
            ;;
        p)
            if [[ $OPTARG =~ ^[0-9]+-[0-9]+$ ]]; then
                p_range="$OPTARG"
            else
                echo "Error: Invalid format for -p. Use start-end (e.g., 377-384)."
                exit 1
            fi
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
    esac
done

# Ensure all ranges are provided
if [[ -z "$a_range" || -z "$b_range" || -z "$p_range" ]]; then
    echo "Error: You must specify all ranges with -a, -b, and -p."
    usage
fi

# Loop through each PDB file
for pdb_file in "$pdb_dir"/*.pdb; do
    # Get the filename without the extension
    pdb_filename=$(basename "$pdb_file" .pdb)

    # Create a directory for the simulation with the same name as the PDB file
    sim_subdir="$pdb_filename"

    # Copy sim_dir scaffold to create simulation directory
    cp -r "$sim_dir" "$sim_subdir"

    # Copy the PDB file to the top_prep directory of the simulation directory
    cp "$pdb_file" "$sim_subdir/top_prep/"

    # Update the tleap file to load the corresponding PDB file
    sed -i "s/mol=loadpdb run_1.pdb/mol=loadpdb $pdb_filename.pdb/" "$sim_subdir/top_prep/tleap.in"

    # Update the bond command in tleap.in
    p1="${p_range%-*}"  # Extract first number from -p
    p2="${p_range#*-}"  # Extract second number from -p
    sed -i "s/bond mol\.[0-9]*\.N mol\.[0-9]*\.C/bond mol.${p1}.N mol.${p2}.C/" "$sim_subdir/top_prep/tleap.in"

    # Update the job name in the submit.sh file
    sed -i "s/#SBATCH --job-name=test_name/#SBATCH --job-name=$pdb_filename/" "$sim_subdir/mdin/submit.sh"

    # Update the atom indices in single.sh command inside submit.sh
    sed -i "s|bash ./single.sh -a [0-9]*-[0-9]* -b [0-9]*-[0-9]* -p [0-9]*-[0-9]*|bash ./single.sh -a ${a_range} -b ${b_range} -p ${p_range}|" "$sim_subdir/mdin/submit.sh"

done
