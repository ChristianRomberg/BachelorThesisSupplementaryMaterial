#!/usr/bin/env python3
import bisect

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import StrMethodFormatter

import paths

preview_factor = 1 # use higher values to skip samples
total_uniref90_gene_families = 87296736  # copied from diamond dbinfo output
aa_counts = pd.read_csv(paths.gene_family_counts, index_col=0)
times_per_split = pd.read_csv(paths.evaluation_times_by_split)

# read total unaligned reads
total_counts_df = pd.read_csv(paths.counts_summary, index_col=0)
total = int(total_counts_df.loc["nucleotide_unaligned"])  # total unaligned reads after nucleotide alignment
unaligned = int(total_counts_df.loc["translated_unaligned"])  # total unaligned reads after translated alignment

total_alignments = aa_counts["count"].sum()
db_alignments_cumsum = np.array(aa_counts["count"]).cumsum()
db_size_fraction = (np.arange(len(db_alignments_cumsum)) + 1) / total_uniref90_gene_families

reads_fraction_aligned_after_first_stage = db_alignments_cumsum / (total_alignments + unaligned)
time_first_stage = db_size_fraction
time_second_stage = (1 - reads_fraction_aligned_after_first_stage) * (1 - db_size_fraction)
time_both_stages = time_first_stage + time_second_stage
quantiles_y = [i * total_alignments / 10 for i in range(1, 10)]
quantiles_x = [bisect.bisect_left(db_alignments_cumsum, quantile) for quantile in quantiles_y]

fig, (plot1, plot2, plot3) = plt.subplots(nrows=3, sharex=True, figsize=(7, 4.5))
# plt.semilogx(cumsum_fraction)
plot1.plot(db_alignments_cumsum[::preview_factor], label="reads covered", c="tab:purple")
for quantile in quantiles_y:
    plot1.axhline(quantile, linestyle=":", color="black", alpha=0.3)
for quantile in quantiles_x:
    plot1.axvline(quantile / preview_factor, linestyle="--", color="black", alpha=0.3)
# plot1.hlines(quantiles_y, -len(cumsum_fraction), 2*len(cumsum_fraction), color="black", alpha=0.3)
# plot1.vlines(quantiles_x, -1, 1, color="black", alpha=0.3)

plot2.plot(time_both_stages[::preview_factor] * 100, label="% time for both stages", c="tab:orange")
plot2.plot(time_first_stage[::preview_factor] * 100, label="% time for first stage", c="tab:green")
plot2.plot(time_second_stage[::preview_factor] * 100, label="% time for second stage", c="tab:blue")

plot3.scatter(times_per_split["split"] / preview_factor, times_per_split["translated alignment"], marker=",", s=2,
              label="time for both stages (s)", c="tab:orange")
plot3.scatter(times_per_split["split"] / preview_factor, times_per_split["translated1"], marker=",", s=2,
              label="time for first stage (s)", c="tab:green")
plot3.scatter(times_per_split["split"] / preview_factor, times_per_split["translated2"], marker=",", s=2,
              label="time for second stage (s)", c="tab:blue")

fig.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
plot3.set_xlabel("Number of gene families")
# plt.ylabel("Percentage of mapped reads")
# plt.title("Gene family vs cumulative abundance")
plot1.legend(loc="upper right")
plot2.legend(loc="upper right")
plot3.legend(loc="upper right")

fig.tight_layout()
# plt.show()
plt.savefig(paths.mathematical_model_vs_reality, dpi=300)

print("Quantiles:")
for quantile, gene_families in zip(quantiles_y, quantiles_x):
    print(f"{quantile / total_alignments:.0%}: {gene_families}")

print("Maximum time saving:")
min_time = np.min(time_both_stages)
print(f"{min_time:.2%}")
print(f"at {np.where(time_both_stages == min_time)[0][0]} gene families")
