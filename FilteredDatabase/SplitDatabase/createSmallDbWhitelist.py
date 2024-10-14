#!/usr/bin/env python3
import bisect

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import StrMethodFormatter

import paths

total_uniref90_gene_families = 87296736  # copied from diamond dbinfo output
aa_counts = pd.read_csv(paths.gene_family_counts, index_col=0)

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

min_time = np.min(time_both_stages)
database_split = np.where(time_both_stages == min_time)[0][0]

print("time taken at min: ", min_time)
print("database split: ", database_split)

small_db_whitelist = set(aa_counts.index[:database_split])

with paths.small_db_genefamilies.open("w") as f:
    f.writelines(s + "\n" for s in sorted(small_db_whitelist))
