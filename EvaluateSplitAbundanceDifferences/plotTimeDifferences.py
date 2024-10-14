#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

import paths

times = pd.read_csv(paths.evaluation_times_by_sample, index_col=[0, 1])

combined_timestamp = "combined nucleotide alignment and post-processing"
slow_indices = times.loc[:, combined_timestamp].isna()
times.loc[slow_indices, combined_timestamp] = (
        times.loc[slow_indices, "nucleotide alignment"] + times.loc[
    slow_indices, "nucleotide alignment post-processing"])

times["other"] = times[
    ["prescreen", "database index", "translated alignment post-processing", "computing gene families",
     "computing pathways"]].sum(axis=1)
times["total"] = times[["other", "custom database creation", "combined nucleotide alignment and post-processing",
                        "translated alignment"]].sum(axis=1)

columns_to_plot = [
    ("custom database creation", "custom database\ncreation"),
    ("combined nucleotide alignment and post-processing", "nucl. alignment\n+ post-processing"),
    ("translated alignment", "transl. alignment"),
    ("other", "other"),
    ("total", "total")
]

boxplot_slow_positions = []
boxplot_slow_values = []
boxplot_fast_positions = []
boxplot_fast_values = []
boxplot_label_positions = []
boxplot_labels = []
current_position = 1
print("Statistics:")

for column, label in columns_to_plot:
    slow_times = times.loc[slow_indices, column].droplevel(1)
    fast_times = times.loc[~slow_indices, column].droplevel(1)

    boxplot_slow_values.append(slow_times)
    boxplot_fast_values.append(fast_times)
    boxplot_slow_positions.append(current_position)
    boxplot_fast_positions.append(current_position + 1)

    boxplot_label_positions.append(current_position + 0.5)
    boxplot_labels.append(label)
    current_position += 3

    print()
    print(column)
    print(f"slow: {slow_times.mean():.1f}s ({slow_times.std():.1f}s)")
    print(f"fast: {fast_times.mean():.1f}s ({fast_times.std():.1f}s)")
    differences = slow_times - fast_times
    print(f"time saved: {differences.mean():.1f}s ({differences.std():.1f}s)")
    print(f"percentage: {(differences / slow_times).mean():.2%} ({(differences / slow_times).std():.2%})")

plt.figure(figsize=(7, 4))
props = dict(color='lightcoral', markeredgecolor='lightcoral', linewidth=1.5)
plt.boxplot(boxplot_slow_values, positions=boxplot_slow_positions, medianprops=dict(color="firebrick"),
            boxprops=props, whiskerprops=props, flierprops=props, capprops=props,
            label="Time (s) for unmodified\nHUMAnN3 version")
props = dict(color='lightskyblue', markeredgecolor='lightskyblue', linewidth=1.5)
plt.boxplot(boxplot_fast_values, positions=boxplot_fast_positions, medianprops=dict(color="dodgerblue"),
            boxprops=props, whiskerprops=props, flierprops=props, capprops=props,
            label="Time (s) for modified\nHUMAnN3 version")
plt.xticks(boxplot_label_positions, boxplot_labels, rotation=15)
plt.ylabel("Time (s)")
plt.legend()
plt.tight_layout()
# plt.show()

plt.savefig("times_pipeline_steps_boxplot.png", dpi=300)
