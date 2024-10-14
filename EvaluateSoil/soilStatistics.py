#!/usr/bin/env python3

import lzma
from pathlib import Path
import pandas as pd
from scipy.stats import pearsonr

import paths


def extract_gene_family_abundance(genefamilies_path):
    def read_file():
        with genefamilies_path.open() as i:
            for line in i:
                name, abundance = line.split("\t")
                if name.endswith("|unclassified"):
                    name = name.split("|")[0]
                    abundance = float(abundance)
                    yield name, abundance

    return (pd.DataFrame(read_file(), columns=["gene_family", "abundance"])
            .set_index("gene_family", drop=True))


small_db_genefamilies = set(paths.small_db_genefamilies.read_text().split())

soil_slow_abundances = extract_gene_family_abundance(paths.soil_slow)
soil_fast_abundances = extract_gene_family_abundance(paths.soil_fast)

soil_slow_abundances.columns = ["slow"]
soil_fast_abundances.columns = ["fast"]

joined = pd.concat([soil_slow_abundances, soil_fast_abundances], axis=1)

soil_slow_abundances = joined.loc[:, "slow"]
soil_fast_abundances = joined.loc[:, "fast"]

small_db_index = joined.index.isin(small_db_genefamilies)


def unmatched_fraction(x, y):
    return y[~x.isna()].isna().sum() / x.count()


def correlation(x, y):
    indices = ~(x.isna() | y.isna())
    return pearsonr(x[indices], y[indices])[0]


print("Total number of gene families in slow")
print(soil_slow_abundances.count())

print("Gene families in slow that were not found in fast")
print("{:.2%}".format(unmatched_fraction(soil_slow_abundances, soil_fast_abundances)))

print("for small db")
print("{:.2%}".format(unmatched_fraction(soil_slow_abundances[small_db_index], soil_fast_abundances[small_db_index])))

print("for large db")
print("{:.2%}".format(unmatched_fraction(soil_slow_abundances[~small_db_index], soil_fast_abundances[~small_db_index])))

print("Total number of gene families in fast")
print(soil_fast_abundances.count())

print("Gene families in fast that were not found in slow")
print("{:.2%}".format(unmatched_fraction(soil_fast_abundances,soil_slow_abundances)))

print("for small db")
print("{:.2%}".format(unmatched_fraction(soil_fast_abundances[small_db_index], soil_slow_abundances[small_db_index])))

print("for large db")
print("{:.2%}".format(unmatched_fraction(soil_fast_abundances[~small_db_index], soil_slow_abundances[~small_db_index])))

print("pearson")
print("{:.4f}".format(correlation(soil_slow_abundances, soil_fast_abundances)))