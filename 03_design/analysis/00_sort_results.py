#!/usr/bin/env python3
"""
Step 00: Parse design run logs and produce results.csv ranked by pLDDT.

Run from 03_design/analysis/. Reads ../design_runs/ by default.
"""

import os
import csv
import re
import argparse


def extract_data_from_log(log_path):
    data = []
    run_id = os.path.basename(os.path.dirname(log_path))
    with open(log_path, 'r') as f:
        for line in f:
            match = re.match(
                r"(\d+) models \[(\d+)\] recycles (\d+) hard (\d+) soft ([\d.]+) "
                r"temp ([\d.]+) loss ([\d.]+) i_con_1 ([\d.]+) i_con_2 ([\d.]+) "
                r"plddt ([\d.]+) ptm ([\d.]+) i_ptm ([\d.]+)", line)
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


def main():
    parser = argparse.ArgumentParser(description='Parse design run logs into results.csv')
    parser.add_argument('--runs_dir', default='../design_runs',
                        help='Design runs directory (default: ../design_runs)')
    parser.add_argument('--output', default='results.csv',
                        help='Output CSV (default: results.csv)')
    args = parser.parse_args()

    all_data = []
    for run_dir in sorted(os.listdir(args.runs_dir)):
        log_path = os.path.join(args.runs_dir, run_dir, "run.log")
        if os.path.isfile(log_path):
            all_data.extend(extract_data_from_log(log_path))

    all_data.sort(key=lambda x: x['plddt'], reverse=True)

    fieldnames = ["run_id", "models", "index", "recycles", "hard", "soft", "temp",
                  "loss", "i_con_1", "i_con_2", "plddt", "ptm", "i_ptm"]
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)

    print(f"Written {len(all_data)} entries to {args.output}")


if __name__ == "__main__":
    main()
