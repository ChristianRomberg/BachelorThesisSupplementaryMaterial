import pandas as pd
import random

import paths

samples_per_study = 10
train_studies = 10
test_studies = 3

metagenomic_data = pd.read_csv("curatedMetagenomicMetadata.csv", low_memory=False, index_col=0)
# filter out samples that we can't access
metagenomic_data = metagenomic_data[~metagenomic_data.NCBI_accession.isna()]

data_filter = metagenomic_data.body_site == "stool"
data_filter &= metagenomic_data.disease == "healthy"
data_filter &= metagenomic_data.age_category == "adult"

# filter out samples with multiple accessions
data_filter &= ~metagenomic_data.NCBI_accession.str.contains(";")

filtered_samples = metagenomic_data[data_filter]

study_sample_counts = filtered_samples.groupby("study_name").count().iloc[:,0]
study_names = [name for name, count in zip(study_sample_counts.index, study_sample_counts) if count >= samples_per_study]

random.seed(42)
studies_for_testing = random.sample(study_names, test_studies)

studies_for_training = [name for name in study_names if name not in studies_for_testing]
studies_for_training = random.sample(studies_for_training, train_studies)

print("Studies for training:")
print(studies_for_training)
print("Accessions for training:")
accessions_for_training = []
for study in studies_for_training:
    study_accessions = list(filtered_samples[filtered_samples.study_name == study].NCBI_accession)
    accessions_for_training.extend(random.sample(study_accessions, samples_per_study))
print(accessions_for_training)
filtered_samples[filtered_samples.NCBI_accession.isin(accessions_for_training)].to_csv(paths.cohort1_metadata)

print("Studies for testing:")
print(studies_for_testing)
print("Accessions for testing:")
accessions_for_testing = []
for study in studies_for_testing:
    study_accessions = list(filtered_samples[filtered_samples.study_name == study].NCBI_accession)
    accessions_for_testing.extend(random.sample(study_accessions, samples_per_study))
print(accessions_for_testing)

filtered_samples[filtered_samples.NCBI_accession.isin(accessions_for_testing)].to_csv(paths.cohort2_metadata)