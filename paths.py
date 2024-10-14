#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).parent
plots = root / "plots"
data = root / "data"

all_genefamilies_abundance = plots / "abundance_all_genefamilies.png"
large_db_abundance = plots / "abundance_large_db.png"
small_db_abundance = plots / "abundance_small_db.png"
pipeline_steps_boxplot = plots / "times_pipeline_steps-boxplot.png"
mathematical_model_vs_reality = plots / "mathematical_model_vs_reality.png"

cohort1_metadata = data / "cohort1.csv"
cohort1_per_sample_results = data / "cohort1_individual_samples"
counts_summary = data / "counts_summary.csv"
gene_family_counts = data / "gene_family_counts.csv"
small_db_genefamilies = data / "small_db_genefamilies.txt"

cohort2_metadata = data / "cohort2.csv"
evaluation_gene_family_abundances = data / "evaluation_gene_family_abundances.csv"
evaluation_times_by_sample = data / "evaluation_times.csv"
evaluation_times_by_split = data / "evaluation_times_by_split.csv"

soil_slow = data / "SRR8881007_genefamilies_unmodified.tsv"
soil_fast = data / "SRR8881007_genefamilies_modified.tsv"
