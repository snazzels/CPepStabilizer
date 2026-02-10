#!/bin/bash

# Default residue and ligand ranges
AE='1-276'
BG='277-376'
LI='377-384'

# Parse command-line arguments
if [ $# -eq 0 ]; then
    echo "Error: No flags provided. Use -a, -b, and -p to specify residue ranges."
    exit 1
fi

while getopts "a:b:p:" opt; do
    case ${opt} in
        a )
            AE="$OPTARG"
            ;;
        b )
            BG="$OPTARG"
            ;;
        p )
            LI="$OPTARG"
            ;;
        \? )
            echo "Usage: $0 [-a A_residues] [-b B_residues] [-p Ligand_residues]"
            exit 1
            ;;
    esac
done

# Path and filenames for your system
prm="../system_wb.prmtop"
nc="../free_run/run.nc"

# Function to process the single system
process_system() {
    local A_res="$AE"
    local B_res="$BG"
    local L_res="$LI"
    local AB_res="${A_res},${B_res}"
    local AL_res="${A_res},${L_res}"
    local LB_res="${B_res},${L_res}"
    local ALB_res="${A_res},${B_res},${L_res}"

    local resids=( $A_res $B_res $L_res $AB_res $AL_res $LB_res $ALB_res )
    local names=( A B L AB AL LB ALB )

    for i in $(seq 0 6); do
        cpptraj=split_${names[$i]}.cpptraj
        echo '' > $cpptraj
        echo "parm $prm" >> $cpptraj
        echo "trajin $nc" >> $cpptraj
        echo "strip !(:${resids[$i]})" >> $cpptraj
        echo "autoimage" >> $cpptraj
        echo "trajout ${names[$i]}.nc" >> $cpptraj
        echo "run" >> $cpptraj
        echo "clear trajin" >> $cpptraj
        echo "parmstrip !(:${resids[$i]})" >> $cpptraj
        echo "parmwrite out ${names[$i]}.prmtop" >> $cpptraj
        cpptraj -i $cpptraj
    done
}

# Call the function to process your system
process_system

