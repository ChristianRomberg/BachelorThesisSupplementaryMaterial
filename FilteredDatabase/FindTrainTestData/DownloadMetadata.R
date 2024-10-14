#!/usr/bin/env R

# R script for acquiring the metadata from curatedMetagenomicData

BiocManager::install("curatedMetagenomicData")
write.csv(curatedMetagenomicData::sampleMetadata, "curatedMetagenomicMetadata.csv")
