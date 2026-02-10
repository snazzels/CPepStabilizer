#!/bin/bash

# Directory where simulations are located
sim_dir="./"

# Find all submit.sh files recursively within the mdin subdirectory of sim_dir
submit_scripts=$(find "$sim_dir" -type f -path "*/mdin/submit.sh")

# Loop through each submit.sh file and execute it
for script in $submit_scripts; do
    # Change directory to the directory containing the script
    script_dir=$(dirname "$script")
    cd "$script_dir" || exit

    # Execute the submit.sh script
    sbatch submit.sh

    # Change back to the original directory
    cd - || exit
done

