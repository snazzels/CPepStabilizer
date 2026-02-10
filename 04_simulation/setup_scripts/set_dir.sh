#!/bin/bash

# Find all run.sh files starting from the current directory
find . -type f -name "run.sh" | while read run_script; do
    # Extract the directory in which run.sh is located
    script_dir=$(dirname "$run_script")

    # Calculate the parent directory of the script directory
    parent_dir=$(dirname "$script_dir")

    # Convert the parent directory path to an absolute path
    absolute_parent_path=$(cd "$parent_dir" && pwd)

    # Prepare the rootdir assignment line with the absolute parent path
    rootdir_line="rootdir=$absolute_parent_path"

    # Check if rootdir is already set in the script
    if grep -q "^rootdir=" "$run_script"; then
        # rootdir is already defined, replace the existing line
        sed -i "/^rootdir=/c\\$rootdir_line" "$run_script"
    else
        # rootdir is not defined, add it after 'set -e'
        sed -i "/^set -e/a\\$rootdir_line" "$run_script"
    fi

    echo "Updated rootdir in $run_script"
done

echo "All scripts have been updated."

