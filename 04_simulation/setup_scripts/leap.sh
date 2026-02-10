#!/bin/bash

# Directory where simulations are located
sim_dir="./"

# Find all run_leap.sh files recursively within the sim_dir
leap_scripts=$(find "$sim_dir" -type f -name "run_leap.sh")

# Loop through each run_leap.sh file and execute it
for script in $leap_scripts; do
    # Change directory to the directory containing the script
    script_dir=$(dirname "$script")
    cd "$script_dir" || exit

    # Execute the run_leap.sh script
    ./run_leap.sh

    # Change back to the original directory
    cd - || exit
done

