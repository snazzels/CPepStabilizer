#!/usr/bin/env python3

import pandas as pd
from scipy.stats import spearmanr, pearsonr, linregress
import matplotlib.pyplot as plt
import os

# Input file
csv_file = "results.csv"

# Load data
df = pd.read_csv(csv_file)

# Define threshold
threshold = 15

# Split data
outliers = df[df["Perplexity"] > threshold]
filtered = df[df["Perplexity"] <= threshold]

print("Outliers (Perplexity > 15):")
print(outliers)
outliers.to_csv("perplexity_outliers.csv", index=False)

# Create figure with 2 rows and 2 columns
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# -------- FULL DATASET --------
x = df["Perplexity"]
y = df["Avg Exact Identity"]

spearman_corr = spearmanr(x, y).correlation
pearson_corr = pearsonr(x, y)[0]
slope, intercept, _, _, _ = linregress(x, y)

# Scatter plot
axes[0, 0].scatter(x, y, alpha=0.6, color="tab:blue", label="Data")
axes[0, 0].plot(x, slope*x + intercept, color="red", linestyle="--", label="Fit")

# Highlight outliers
axes[0, 0].scatter(outliers["Perplexity"], outliers["Avg Exact Identity"],
                   color="red", edgecolor="black", s=100, label="Outliers")

if "PDB_ID" in outliers.columns:
    for _, row in outliers.iterrows():
        axes[0, 0].annotate(row["PDB_ID"], 
                            (row["Perplexity"], row["Avg Exact Identity"]),
                            textcoords="offset points", xytext=(5,5), fontsize=8)

# Threshold line
axes[0, 0].axvline(x=threshold, color="gray", linestyle=":", label=f"Threshold={threshold}")

axes[0, 0].set_xlabel("Perplexity")
axes[0, 0].set_ylabel("Avg Exact Identity")
axes[0, 0].set_title("Perplexity vs Avg Exact Identity (All samples)")
axes[0, 0].legend(title=f"ρ={spearman_corr:.2f}, r={pearson_corr:.2f}")
axes[0, 0].grid(True)

# Boxplot
axes[0, 1].boxplot(y, patch_artist=True)
axes[0, 1].set_ylabel("Avg Exact Identity")
axes[0, 1].set_title("Sequence Identity Distribution (All samples)")
axes[0, 1].set_xticklabels(["output"])
axes[0, 1].grid(True)

# -------- FILTERED DATASET --------
x_f = filtered["Perplexity"]
y_f = filtered["Avg Exact Identity"]

if len(filtered) > 1:
    spearman_corr_f = spearmanr(x_f, y_f).correlation
    pearson_corr_f = pearsonr(x_f, y_f)[0]
    slope_f, intercept_f, _, _, _ = linregress(x_f, y_f)
else:
    spearman_corr_f = pearson_corr_f = float("nan")
    slope_f, intercept_f = 0, y_f.mean() if len(y_f) else 0

# Scatter plot
axes[1, 0].scatter(x_f, y_f, alpha=0.6, color="tab:green", label="Data")
if len(filtered) > 1:
    axes[1, 0].plot(x_f, slope_f*x_f + intercept_f, color="red", linestyle="--", label="Fit")

axes[1, 0].set_xlabel("Perplexity")
axes[1, 0].set_ylabel("Avg Exact Identity")
axes[1, 0].set_title("Perplexity vs Avg Exact Identity (Perplexity ≤ 15)")
axes[1, 0].legend(title=f"ρ={spearman_corr_f:.2f}, r={pearson_corr_f:.2f}")
axes[1, 0].grid(True)

# Boxplot
axes[1, 1].boxplot(y_f, patch_artist=True)
axes[1, 1].set_ylabel("Avg Exact Identity")
axes[1, 1].set_title("Sequence Identity Distribution (Perplexity ≤ 15)")
axes[1, 1].set_xticklabels(["output"])
axes[1, 1].grid(True)

plt.tight_layout()
plt.show()

