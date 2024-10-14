#!/usr/bin/env python3
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, pearsonr

import paths

abundances = pd.read_csv(paths.evaluation_gene_family_abundances, index_col=0)
abundances.rename(columns=eval, inplace=True)

small_db_genefamilies = set(paths.small_db_genefamilies.read_text().split())

accessions, types = zip(*abundances.columns)


def plot(x, y, outfile_name=None):
    x = np.log2(x.fillna(0) + 1)
    y = np.log2(y.fillna(0) + 1)

    xy = list(zip(x, y))
    xy = np.array(random.sample(xy, 100000))
    x, y = xy[:, 0], xy[:, 1]

    c = gaussian_kde(xy.transpose())(xy.transpose())

    # plt.hist2d(x, y, bins=1000)
    plt.figure(figsize=(2.5, 2.5))
    plt.scatter(x, y, c=c, s=0.1, marker=',')

    plt.xlabel("log2 abundance")
    plt.ylabel("log2 abundance")
    plt.xlim(-.2, 10)
    plt.ylim(-.2, 10)

    if outfile_name is not None:
        plt.savefig(outfile_name, bbox_inches='tight', dpi=300)
    else:
        plt.show()


all_genefamilies_x = []
all_genefamilies_y = []

small_db_genefamilies_x = []
small_db_genefamilies_y = []

large_db_genefamilies_x = []
large_db_genefamilies_y = []

for accession in accessions:
    print(accession)
    if (accession, "fast") in abundances.columns and (accession, "slow") in abundances.columns:
        sub_df = abundances[[(accession, "slow"), (accession, "fast")]]
        sub_df = sub_df.dropna(axis=0, how="all")

        x = sub_df.iloc[:, 0]
        y = sub_df.iloc[:, 1]
        all_genefamilies_x.append(x)
        all_genefamilies_y.append(y)

        small_indices = sub_df.index.isin(small_db_genefamilies)

        x_small = x[small_indices]
        y_small = y[small_indices]
        small_db_genefamilies_x.append(x_small)
        small_db_genefamilies_y.append(y_small)
        # plot(x_small, y_small)

        x_large = x[~small_indices]
        y_large = y[~small_indices]
        large_db_genefamilies_x.append(x_large)
        large_db_genefamilies_y.append(y_large)
        # plot(x_large, y_large)
        # break


def unmatched_fraction(x, y):
    return y[~x.isna()].isna().sum() / x.count()


def correlation(x, y):
    indices = ~(x.isna() | y.isna())
    return pearsonr(x[indices], y[indices])[0]


print("Statistics:")
print("Percentage of slow found gene families that didn't have a corresponding gene family in split search")
percentages = [unmatched_fraction(x, y) for x, y in zip(all_genefamilies_x, all_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("for gene families in small db")
percentages = [unmatched_fraction(x, y) for x, y in zip(small_db_genefamilies_x, small_db_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("for gene families in large db")
percentages = [unmatched_fraction(x, y) for x, y in zip(large_db_genefamilies_x, large_db_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("Percentage of split found gene families that didn't have a corresponding gene family in slow search")
percentages = [unmatched_fraction(y, x) for x, y in zip(all_genefamilies_x, all_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("for gene families in small db")
percentages = [unmatched_fraction(y, x) for x, y in zip(small_db_genefamilies_x, small_db_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("for gene families in large db")
percentages = [unmatched_fraction(y, x) for x, y in zip(large_db_genefamilies_x, large_db_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("Pearson correllation coefficient")
values = [correlation(x, y) for x, y in zip(all_genefamilies_x, all_genefamilies_y)]
print(f"Average: {np.mean(values):.4f}; std: {np.std(values):.4f}")

print("for gene families in small db")
values = [correlation(x, y) for x, y in zip(small_db_genefamilies_x, small_db_genefamilies_y)]
print(f"Average: {np.mean(values):.4f}; std: {np.std(values):.4f}")

print("for gene families in large db")
values = [correlation(x, y) for x, y in zip(large_db_genefamilies_x, large_db_genefamilies_y)]
print(f"Average: {np.mean(values):.4f}; std: {np.std(values):.4f}")

print("Percentage of slow found genefamilies that are in large db")
percentages = [(~large.isna()).sum() / (~all_.isna()).sum() for all_, large in
               zip(all_genefamilies_x, large_db_genefamilies_x)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("Percentage of slow found genefamilies that are in small db")
percentages = [(~small.isna()).sum() / (~all_.isna()).sum() for all_, small in
               zip(all_genefamilies_x, small_db_genefamilies_x)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("Percentage of split found genefamilies that are in large db")
percentages = [(~large.isna()).sum() / (~all_.isna()).sum() for all_, large in
               zip(all_genefamilies_y, large_db_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("Percentage of split found genefamilies that are in small db")
percentages = [(~small.isna()).sum() / (~all_.isna()).sum() for all_, small in
               zip(all_genefamilies_y, small_db_genefamilies_y)]
print(f"Average: {np.mean(percentages):.2%}; std: {np.std(percentages):.2%}")

print("Number of gene families present in slow search")
values = [(~x.isna()).sum() for x in all_genefamilies_x]
print(f"Average: {np.mean(values):.4f}; std: {np.std(values):.4f}")

print("Number of gene families present in split search")
values = [(~y.isna()).sum() for y in all_genefamilies_y]
print(f"Average: {np.mean(values):.4f}; std: {np.std(values):.4f}")

print("Starting plots")
plot(pd.concat(all_genefamilies_x, ignore_index=True), pd.concat(all_genefamilies_y, ignore_index=True),
     "abundance_all_genefamilies.png")
plot(pd.concat(small_db_genefamilies_x, ignore_index=True), pd.concat(small_db_genefamilies_y, ignore_index=True),
     "abundance_small_db_genefamilies.png")
plot(pd.concat(large_db_genefamilies_x, ignore_index=True), pd.concat(large_db_genefamilies_y, ignore_index=True),
     "abundance_large_db_genefamilies.png")
