#!/usr/bin/env python3
"""Collect MM/GBSA results from all run_sim directories, merge with design
metrics, compute cooperativity, and save a summary CSV + plots.

Output (in run_sim/):
  mmgbsa_summary.csv  — merged design + MMGBSA metrics, ranked by dG_AB_L
  mmgbsa_plots.png    — correlation plots (BSA, pLDDT, cooperativity)

Calculation meanings:
  AL   — peptide (L) binding to chain A alone
  BL   — peptide (L) binding to chain B alone
  AB_L — peptide (L) binding to the A+B complex  (= dG_bind of peptide)
  AB   — chain A binding to chain B (dimer stability without peptide)
  AL_B — chain A+L complex binding to chain B   (dimer stability with peptide)

Cooperativity:
  coop = dG_AL_B - dG_AB
  Negative → peptide stabilises the dimer; positive → peptide disrupts it.
  (Equivalent by thermodynamic cycle to dG_AB_L - dG_AL.)
"""

import os
import re
import sys
import pathlib
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GBSA_PATTERN = re.compile(
    r'DELTA TOTAL\s+(?P<avg>-?\d+\.\d+)\s+(?P<std>\d+\.\d+)\s+(?P<sem>\d+\.\d+)'
)

CALC_TYPES = [
    ("AL",   "AL_out/AL_output.dat"),
    ("BL",   "BL_out/BL_output.dat"),
    ("AB_L", "AB_L_out/AB_L_output.dat"),
    ("AB",   "AB_out/AB_output.dat"),
    ("AL_B", "AL_B_out/AL_B_output.dat"),
]


def parse_output(path):
    """Return (average, std, sem) from an MMPBSA output file, or None."""
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        content = f.read()
    m = GBSA_PATTERN.search(content)
    if m is None:
        return None
    return float(m.group("avg")), float(m.group("std")), float(m.group("sem"))


def collect(run_sim_dir):
    run_sim_dir = pathlib.Path(run_sim_dir)
    rows = []

    for run_dir in sorted(run_sim_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        mmgbsa_dir = run_dir / "mmgbsa"
        if not mmgbsa_dir.is_dir():
            continue

        row = {"run_name": run_dir.name}
        any_found = False

        for label, rel_path in CALC_TYPES:
            result = parse_output(mmgbsa_dir / rel_path)
            if result is not None:
                avg, std, sem = result
                any_found = True
            else:
                avg, std, sem = None, None, None
            row[f"dG_{label}"]  = avg
            row[f"std_{label}"] = std
            row[f"sem_{label}"] = sem

        if any_found:
            rows.append(row)
        else:
            print(f"  [skip] {run_dir.name} — no completed MMGBSA output found", file=sys.stderr)

    return pd.DataFrame(rows)


def merge_design_results(df_mmgbsa, script_dir):
    """Merge MMGBSA data with design-phase results.csv."""
    results_path = script_dir.parent / "03_design" / "analysis" / "results.csv"
    if not results_path.is_file():
        print(f"  [info] Design results not found at {results_path}, skipping merge.", file=sys.stderr)
        return df_mmgbsa

    design = pd.read_csv(results_path)
    if "pdb_filename" not in design.columns:
        print("  [info] results.csv has no pdb_filename column, skipping merge.", file=sys.stderr)
        return df_mmgbsa

    # Build merge key: strip .pdb suffix to match run_name
    design["run_name"] = design["pdb_filename"].str.replace(r"\.pdb$", "", regex=True)

    keep = [c for c in ["run_name", "plddt", "ptm", "i_ptm", "bsa_a", "bsa_b",
                         "sasa_peptide", "Avg Exact Identity", "Avg Adjusted Identity"]
            if c in design.columns]
    merged = df_mmgbsa.merge(design[keep], on="run_name", how="left")

    if "bsa_a" in merged.columns and "bsa_b" in merged.columns:
        merged["bsa_total"] = merged["bsa_a"] + merged["bsa_b"]
        merged["bsa_min"]   = merged[["bsa_a", "bsa_b"]].min(axis=1)

    n_matched = merged["plddt"].notna().sum()
    print(f"  Merged with design results: {n_matched}/{len(merged)} runs matched.")
    return merged


def compute_cooperativity(df):
    """Add cooperativity column: coop = dG_AL_B - dG_AB."""
    if "dG_AL_B" in df.columns and "dG_AB" in df.columns:
        df["coop"] = df["dG_AL_B"] - df["dG_AB"]
    return df


def make_plots(df, output):
    has_design = "plddt" in df.columns
    has_bsa    = "bsa_min" in df.columns
    has_coop   = "coop" in df.columns

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("MM/GBSA Summary", fontsize=14)

    def scatter(ax, x, y, c=None, xlabel="", ylabel="", title="", cmap="viridis", clabel=""):
        sub = df[[x, y] + ([c] if c else [])].dropna()
        if sub.empty:
            ax.set_visible(False)
            return
        sc = ax.scatter(sub[x], sub[y], c=sub[c] if c else "steelblue",
                        cmap=cmap if c else None, alpha=0.8, s=60, edgecolors="white", linewidths=0.3)
        if c:
            cb = fig.colorbar(sc, ax=ax, pad=0.02)
            cb.set_label(clabel or c, fontsize=8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

    # Row 1: key correlations
    if has_bsa:
        scatter(axes[0, 0], "bsa_min", "dG_AB_L",
                c="plddt" if has_design else None,
                xlabel="min(BSA_A, BSA_B) (Å²)", ylabel="dG_AB_L (kcal/mol)",
                title="Min BSA vs Peptide binding energy",
                clabel="pLDDT")
    else:
        axes[0, 0].set_visible(False)

    if has_design:
        scatter(axes[0, 1], "plddt", "dG_AB_L",
                c="bsa_min" if has_bsa else None,
                xlabel="pLDDT", ylabel="dG_AB_L (kcal/mol)",
                title="pLDDT vs Peptide binding energy",
                cmap="plasma", clabel="min BSA (Å²)")
    else:
        axes[0, 1].set_visible(False)

    scatter(axes[0, 2], "dG_AL", "dG_BL",
            xlabel="dG_AL (kcal/mol)", ylabel="dG_BL (kcal/mol)",
            title="Chain A vs Chain B binding")
    # Add diagonal reference line
    sub_ab = df[["dG_AL", "dG_BL"]].dropna()
    if not sub_ab.empty:
        lim = [min(sub_ab.min()), max(sub_ab.max())]
        axes[0, 2].plot(lim, lim, "k--", linewidth=0.8, alpha=0.5)

    # Row 2: cooperativity + dimer stability
    if has_coop:
        sub_c = df["coop"].dropna()
        if not sub_c.empty:
            axes[1, 0].hist(sub_c, bins=max(5, len(sub_c) // 2),
                            color="mediumpurple", edgecolor="white")
            axes[1, 0].axvline(0, color="black", linewidth=0.8, linestyle="--")
            axes[1, 0].set_xlabel("Cooperativity (kcal/mol)")
            axes[1, 0].set_ylabel("Count")
            axes[1, 0].set_title("Cooperativity: dG_AL_B − dG_AB\n(neg = peptide stabilises dimer)")
        else:
            axes[1, 0].set_visible(False)

        scatter(axes[1, 1], "dG_AB_L", "coop",
                xlabel="dG_AB_L (kcal/mol)", ylabel="Cooperativity (kcal/mol)",
                title="Binding energy vs Cooperativity")
        if not df[["dG_AB_L", "coop"]].dropna().empty:
            axes[1, 1].axhline(0, color="black", linewidth=0.8, linestyle="--")
    else:
        axes[1, 0].set_visible(False)
        axes[1, 1].set_visible(False)

    # dG_AB_L distribution
    sub_abl = df["dG_AB_L"].dropna()
    if not sub_abl.empty:
        axes[1, 2].hist(sub_abl, bins=max(5, len(sub_abl) // 2),
                        color="steelblue", edgecolor="white")
        axes[1, 2].set_xlabel("dG_AB_L (kcal/mol)")
        axes[1, 2].set_ylabel("Count")
        axes[1, 2].set_title("Peptide binding energy distribution")
    else:
        axes[1, 2].set_visible(False)

    plt.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output}")


def main():
    script_dir = pathlib.Path(__file__).parent
    run_sim_dir = script_dir / "run_sim"

    if not run_sim_dir.is_dir():
        print(f"ERROR: {run_sim_dir} not found", file=sys.stderr)
        sys.exit(1)

    print(f"Collecting MM/GBSA results from: {run_sim_dir}")
    df = collect(run_sim_dir)

    if df.empty:
        print("No completed MM/GBSA results found.", file=sys.stderr)
        sys.exit(1)

    df = merge_design_results(df, script_dir)
    df = compute_cooperativity(df)

    # Sort by AB_L binding energy (most negative = strongest binder)
    df = df.sort_values("dG_AB_L", ascending=True, na_position="last").reset_index(drop=True)
    df.index += 1

    out_csv = run_sim_dir / "mmgbsa_summary.csv"
    df.to_csv(out_csv, index_label="rank")
    print(f"Saved: {out_csv}")

    make_plots(df, run_sim_dir / "mmgbsa_plots.png")
    print()

    # Print ranked table
    display_cols = ["run_name", "dG_AB_L", "std_AB_L", "coop",
                    "dG_AL", "dG_BL", "dG_AB"]
    display_cols = [c for c in display_cols if c in df.columns]
    if "plddt" in df.columns:
        display_cols.insert(1, "plddt")
    if "bsa_min" in df.columns:
        display_cols.insert(2, "bsa_min")

    print(df[display_cols].to_string(
        float_format=lambda x: f"{x:7.2f}" if pd.notna(x) else "    N/A"
    ))


if __name__ == "__main__":
    main()
