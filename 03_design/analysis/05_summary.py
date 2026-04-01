#!/usr/bin/env python3
"""
Step 08: Print ranked table of top designs and save summary plots.

Reads results.csv and produces:
  - Terminal output: ranked table of top 10 designs
  - summary_plots.png: 2x4 panel with histograms (row 1) and scatter plots (row 2)
"""

import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_results(csv_path="results.csv"):
    df = pd.read_csv(csv_path)
    if "bsa_a" in df.columns and "bsa_b" in df.columns:
        df["bsa_total"] = df["bsa_a"] + df["bsa_b"]
        df["bsa_min"] = df[["bsa_a", "bsa_b"]].min(axis=1)
    return df


def print_top_table(df, n=10):
    required = {"pdb_filename", "plddt"}
    if not required.issubset(df.columns):
        print("results.csv is missing required columns (pdb_filename, plddt).", file=sys.stderr)
        return

    optional = ["ptm", "i_ptm", "bsa_total", "bsa_min", "Avg Exact Identity"]
    display_cols = ["pdb_filename", "plddt"] + [c for c in optional if c in df.columns]

    top = df.sort_values("plddt", ascending=False).head(n)[display_cols].reset_index(drop=True)
    top.index += 1

    col_headers = {
        "pdb_filename": "pdb_filename",
        "plddt": "pLDDT",
        "ptm": "ptm",
        "i_ptm": "i_ptm",
        "bsa_total": "BSA_total",
        "bsa_min": "BSA_min",
        "Avg Exact Identity": "Seq_ID",
    }
    top = top.rename(columns=col_headers)

    # Format floats
    for col in ["pLDDT", "ptm", "i_ptm", "Seq_ID"]:
        if col in top.columns:
            top[col] = top[col].map(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    for col in ["BSA_total", "BSA_min"]:
        if col in top.columns:
            top[col] = top[col].map(lambda x: f"{x:.1f}" if pd.notna(x) else "")

    print("\nTop designs ranked by pLDDT:")
    print(top.to_string())
    print()


def make_summary_plots(df, output="summary_plots.png"):
    has_bsa = "bsa_total" in df.columns and "bsa_min" in df.columns
    has_seqid = "Avg Exact Identity" in df.columns

    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    fig.suptitle("Design Summary", fontsize=14)

    # Row 1: histograms
    axes[0, 0].hist(df["plddt"].dropna(), bins=20, color="steelblue", edgecolor="white")
    axes[0, 0].set_xlabel("pLDDT")
    axes[0, 0].set_ylabel("Count")
    axes[0, 0].set_title("pLDDT distribution")

    if has_bsa:
        axes[0, 1].hist(df["bsa_total"].dropna(), bins=20, color="darkorange", edgecolor="white")
        axes[0, 1].set_xlabel("BSA total (Å²)")
        axes[0, 1].set_title("Total BSA distribution")

        axes[0, 2].hist(df["bsa_min"].dropna(), bins=20, color="mediumpurple", edgecolor="white")
        axes[0, 2].set_xlabel("min(BSA_A, BSA_B) (Å²)")
        axes[0, 2].set_title("Min BSA distribution")
    else:
        axes[0, 1].set_visible(False)
        axes[0, 2].set_visible(False)

    if has_seqid:
        axes[0, 3].hist(df["Avg Exact Identity"].dropna(), bins=20, color="seagreen", edgecolor="white")
        axes[0, 3].set_xlabel("Avg Exact Identity")
        axes[0, 3].set_title("Sequence identity distribution")
    else:
        axes[0, 3].set_visible(False)

    # Row 2: scatter plots
    if has_bsa:
        sub = df[["plddt", "bsa_total"]].dropna()
        axes[1, 0].scatter(sub["plddt"], sub["bsa_total"], alpha=0.5, s=15, color="steelblue")
        axes[1, 0].set_xlabel("pLDDT")
        axes[1, 0].set_ylabel("BSA total (Å²)")
        axes[1, 0].set_title("pLDDT vs BSA total")

        sub = df[["plddt", "bsa_min"]].dropna()
        axes[1, 1].scatter(sub["plddt"], sub["bsa_min"], alpha=0.5, s=15, color="mediumpurple")
        axes[1, 1].set_xlabel("pLDDT")
        axes[1, 1].set_ylabel("min(BSA_A, BSA_B) (Å²)")
        axes[1, 1].set_title("pLDDT vs min BSA")
    else:
        axes[1, 0].set_visible(False)
        axes[1, 1].set_visible(False)

    if has_seqid:
        sub = df[["plddt", "Avg Exact Identity"]].dropna()
        axes[1, 2].scatter(sub["plddt"], sub["Avg Exact Identity"], alpha=0.5, s=15, color="darkorange")
        axes[1, 2].set_xlabel("pLDDT")
        axes[1, 2].set_ylabel("Avg Exact Identity")
        axes[1, 2].set_title("pLDDT vs Seq identity")
    else:
        axes[1, 2].set_visible(False)

    if has_bsa and has_seqid:
        sub = df[["bsa_min", "Avg Exact Identity"]].dropna()
        axes[1, 3].scatter(sub["bsa_min"], sub["Avg Exact Identity"], alpha=0.5, s=15, color="seagreen")
        axes[1, 3].set_xlabel("min(BSA_A, BSA_B) (Å²)")
        axes[1, 3].set_ylabel("Avg Exact Identity")
        axes[1, 3].set_title("min BSA vs Seq identity")
    else:
        axes[1, 3].set_visible(False)

    plt.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main():
    df = load_results()
    print_top_table(df)
    make_summary_plots(df)


if __name__ == "__main__":
    main()
