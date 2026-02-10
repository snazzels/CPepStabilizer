import os
import csv
import re

def extract_data_from_log(log_path):
    data = []
    run_id = os.path.basename(os.path.dirname(log_path))  # Extract run directory name
    
    with open(log_path, 'r') as f:
        for line in f:
            match = re.match(r"(\d+) models \[(\d+)\] recycles (\d+) hard (\d+) soft ([\d.]+) temp ([\d.]+) loss ([\d.]+) i_con_1 ([\d.]+) i_con_2 ([\d.]+) plddt ([\d.]+) ptm ([\d.]+) i_ptm ([\d.]+)", line)
            if match:
                data.append({
                    "run_id": run_id,
                    "models": int(match.group(1)),
                    "index": int(match.group(2)),
                    "recycles": int(match.group(3)),
                    "hard": int(match.group(4)),
                    "soft": float(match.group(5)),
                    "temp": float(match.group(6)),
                    "loss": float(match.group(7)),
                    "i_con_1": float(match.group(8)),
                    "i_con_2": float(match.group(9)),
                    "plddt": float(match.group(10)),
                    "ptm": float(match.group(11)),
                    "i_ptm": float(match.group(12)),
                })
    return data

def process_logs(root_dir, output_csv):
    all_data = []
    
    for run_dir in sorted(os.listdir(root_dir)):
        log_path = os.path.join(root_dir, run_dir, "run.log")
        if os.path.isfile(log_path):
            all_data.extend(extract_data_from_log(log_path))
    
    all_data.sort(key=lambda x: x['plddt'], reverse=True)  # Sort by plddt descending
    
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ["run_id", "models", "index", "recycles", "hard", "soft", "temp", "loss", "i_con_1", "i_con_2", "plddt", "ptm", "i_ptm"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)

# Example usage:
process_logs("design_runs/", "sorted_results.csv")

