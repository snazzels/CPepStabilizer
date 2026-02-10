#!/usr/bin/env python3


import os
import shutil
import pandas as pd

def process_design_results(input_csv, output_filtered_csv, best_top_dir, design_runs_dir):
    """
    Process design results by filtering and copying best PDB files.
    
    Parameters:
    - input_csv: Path to the input CSV file with design results
    - output_filtered_csv: Path to save the filtered and updated CSV
    - best_top_dir: Directory to copy the best PDB files
    - design_runs_dir: Root directory containing design run directories
    """
    # Create output directory if it doesn't exist
    os.makedirs(best_top_dir, exist_ok=True)

    # Load the CSV file
    df = pd.read_csv(input_csv)

    # Filter the dataframe based on specific criteria
    filtered_df = df[(df['hard'] == 1) & (df['plddt'] >= 0.8)]

    # List to store new filenames
    filenames = []

    # Loop through the filtered dataframe
    for _, row in filtered_df.iterrows():
        run_id = row['run_id']  # Get the run_id for the directory
        model_number = row['models'] - 1  # Get the model number (corresponding to model_step_*)
        
        trajectory_dir = os.path.join(design_runs_dir, run_id, 'trajectory')  # Full path to the trajectory directory
        
        # Find the pdb file corresponding to the model number (e.g., model_step_1.pdb)
        pdb_file = os.path.join(trajectory_dir, f'model_step_{model_number}.pdb')

        # Check if the file exists
        if os.path.isfile(pdb_file):
            new_filename = f'{run_id}_model_step_{model_number}.pdb'
            shutil.copy(pdb_file, os.path.join(best_top_dir, new_filename))
            filenames.append(new_filename)
            print(f'Copied {pdb_file} to {best_top_dir}')
        else:
            filenames.append(None)
            print(f'File {pdb_file} not found for run_id {run_id} and model_step_{model_number}')

    # Add the new column to the dataframe
    filtered_df['pdb_filename'] = filenames

    # Save the filtered and updated dataframe
    filtered_df.to_csv(output_filtered_csv, index=False)

    # List files in best_top directory
    print(f'Files in {best_top_dir} directory: {os.listdir(best_top_dir)}')
    print(f'Updated CSV saved as {output_filtered_csv}')

    return filtered_df

def main():
    # Configuration parameters
    input_csv = 'sorted_results.csv'
    output_filtered_csv = 'filtered_best.csv'
    best_top_dir = 'best_top'
    design_runs_dir = '../design_runs'

    # Process the design results
    result_df = process_design_results(
        input_csv, 
        output_filtered_csv, 
        best_top_dir, 
        design_runs_dir
    )

    # Optional: Print the first few rows of the filtered dataframe
    print("\nFiltered Results:")
    print(result_df)

if __name__ == "__main__":
    main()
