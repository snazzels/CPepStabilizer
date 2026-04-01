#!/usr/bin/env python3

import os
import argparse
import subprocess
import tempfile
import pandas as pd

def add_ter_lines(pdb_content):
    """Add TER lines between residues with non-sequential IDs"""
    lines = pdb_content.splitlines(True)
    new_lines = []
    last_residue_id = None

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Extract residue ID from the line
            current_residue_id = int(line[22:26].strip())

            if last_residue_id is not None and current_residue_id != last_residue_id:
                # If there is a jump in residue IDs, add a TER line
                if current_residue_id != last_residue_id + 1:
                    new_lines.append("TER\n")

            last_residue_id = current_residue_id
            new_lines.append(line)
        else:
            # Add lines that are not ATOM or HETATM as they are
            new_lines.append(line)

    # Ensure there is a TER line at the end if the last line was an ATOM/HETATM line
    if new_lines and not new_lines[-1].strip().startswith("TER"):
        new_lines.append("TER\n")

    return ''.join(new_lines)

def assign_chain_identifiers(pdb_content):
    """Assign chain identifiers to the PDB content"""
    lines = pdb_content.splitlines(True)
    new_lines = []
    current_chain = 'A'  # Start with chain A
    chain_switched = False

    for line in lines:
        if line.startswith("ATOM"):
            # Insert the chain identifier directly into the correct column
            line = line[:21] + current_chain + line[22:]
            new_lines.append(line)
        elif line.startswith("TER"):
            # Write the TER line with the current chain identifier and switch chain
            ter_line = line[:21] + current_chain + line[22:]
            new_lines.append(ter_line)
            if not chain_switched:
                current_chain = 'B'
                chain_switched = True
            else:
                current_chain = chr(ord(current_chain) + 1)  # Increment to the next chain identifier
        else:
            # Write all other lines unchanged
            new_lines.append(line)

    return ''.join(new_lines)

def run_pdb4amber(input_file, output_file):
    """Run pdb4amber on the input file to clean it"""
    cmd = f"pdb4amber -i {input_file} -o {output_file}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"pdb4amber processed: {output_file}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to run pdb4amber on {input_file}: {e}")
        return False

def process_pdb_file(input_file, output_dir, use_pdb4amber=True):
    """Process a single PDB file with both ter lines and chain IDs"""
    base_name = os.path.basename(input_file)
    output_file = os.path.join(output_dir, base_name)

    with open(input_file, 'r') as file:
        pdb_content = file.read()

    # Step 1: Add TER lines
    pdb_content_with_ter = add_ter_lines(pdb_content)

    # Step 2: Clean with pdb4amber if requested (uses system temp dir)
    if use_pdb4amber:
        with tempfile.NamedTemporaryFile(suffix="_ter.pdb", mode='w', delete=False) as tf_ter:
            tf_ter.write(pdb_content_with_ter)
            temp_ter_path = tf_ter.name
        with tempfile.NamedTemporaryFile(suffix="_cleaned.pdb", delete=False) as tf_clean:
            temp_clean_path = tf_clean.name
        try:
            success = run_pdb4amber(temp_ter_path, temp_clean_path)
            if success:
                with open(temp_clean_path, 'r') as file:
                    pdb_content = file.read()
            else:
                pdb_content = pdb_content_with_ter
        finally:
            os.unlink(temp_ter_path)
            if os.path.exists(temp_clean_path):
                os.unlink(temp_clean_path)
    else:
        pdb_content = pdb_content_with_ter

    # Step 3: Assign chain identifiers
    final_pdb_content = assign_chain_identifiers(pdb_content)

    with open(output_file, 'w') as file:
        file.write(final_pdb_content)

    print(f"Processed {input_file} and saved to {output_file}")
    return output_file

def main():
    parser = argparse.ArgumentParser(description="Process PDB files to add TER lines and assign chain identifiers.")
    parser.add_argument("input", nargs="?", default=None,
                        help="PDB file or directory to process. Omit to read file list from results.csv.")
    parser.add_argument("-o", "--output_dir", default="pdbs", help="The directory to save processed PDB files")
    parser.add_argument("--skip_pdb4amber", action="store_true", help="Skip the pdb4amber cleaning step")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.input is None:
        # Default: process only PDBs listed in results.csv
        df = pd.read_csv("results.csv")
        for pdb_filename in df["pdb_filename"].dropna():
            input_file = os.path.join("pdbs", pdb_filename)
            if os.path.isfile(input_file):
                process_pdb_file(input_file, args.output_dir, not args.skip_pdb4amber)
            else:
                print(f"File not found, skipping: {input_file}")
    elif os.path.isdir(args.input):
        for filename in os.listdir(args.input):
            if filename.endswith(".pdb"):
                process_pdb_file(os.path.join(args.input, filename), args.output_dir, not args.skip_pdb4amber)
    elif os.path.isfile(args.input) and args.input.endswith(".pdb"):
        process_pdb_file(args.input, args.output_dir, not args.skip_pdb4amber)
    else:
        print(f"Error: {args.input} is not a valid PDB file or directory")

if __name__ == "__main__":
    main()
