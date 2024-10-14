#!/usr/bin/env python3

import lzma
from pathlib import Path
import pandas as pd

import paths

slow1_path = Path(__file__).parent / "humann_results_slow"
fast_path = Path(__file__).parent / "humann_results_fast"

input_paths = {
    "slow": slow1_path,
    "fast": fast_path,
}

def extract_gene_family_abundance(genefamilies_path):
    def read_file():
        with lzma.open(genefamilies_path, "rt") as i:
            for line in i:
                name, abundance = line.split("\t")
                if name.endswith("|unclassified"):
                    name = name.split("|")[0]
                    abundance = float(abundance)
                    yield name, abundance

    return (pd.DataFrame(read_file(), columns=["gene_family", "abundance"])
            .set_index("gene_family", drop=True))


def extract_all_gene_family_abundances():
    gene_families_dfs = []
    for mode, path in input_paths.items():
        for accession_path in path.iterdir():
            df: pd.DataFrame = extract_gene_family_abundance(accession_path / "genefamilies.tsv.xz")
            df = df.rename(columns={"abundance": (accession_path.name, mode)})
            gene_families_dfs.append(df)
    return pd.concat(gene_families_dfs, axis=1)


complete_df = extract_all_gene_family_abundances()
complete_df.to_csv(paths.evaluation_gene_family_abundances)
