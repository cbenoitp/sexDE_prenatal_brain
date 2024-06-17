#!/bin/bash

# go to working directory

# from the E-MTAB-4840.sdrf.txt file all unique ENA run accession number
cut -f 36 data/0-input/E-MTAB-4840.sdrf.txt | sort | uniq | tail -n +2 > data/0-input/PRJEB14954_samples_id.txt
wc -l data/0-input/PRJEB14954_samples_id.txt
## 650 uniq id have been found

# from counts.txt.gz file extract all lines corresponding to the previously selected samples ## very very long command run on cluster
zfgrep -fdata/0-input/PRJEB14954_samples_id.txt data/0-input/counts.txt.gz > data/1-preprocessing/counts_PRJEB14594_without_header.txt
zcat < data/0-input/counts.txt.gz | head -n1 > data/1-preprocessing/header.txt
cat data/1-preprocessing/header.txt data/1-preprocessing/counts_PRJEB14594_without_header.txt > data/1-preprocessing/counts_PRJEB14594.txt


# from the counts extracted, get the id of the samples with counts
cut -f1 data/1-preprocessing/counts_PRJEB14594.txt | tail -n+2 > data/1-preprocessing/PRJEB14954_samples_id_with_counts.txt
wc -l data/1-preprocessing/PRJEB14954_samples_id_with_counts.txt
## 558 samples with counts

# Extract only the info of interest for the samples with counts
## dev_stage, brain_region, karyotype, sample_name and individual
grep -f data/1-preprocessing/PRJEB14954_samples_id_with_counts.txt data/0-input/E-MTAB-4840.sdrf.txt | cut -f 1,5,6,7,9,11,36 | sort -k4 | uniq > data/1-preprocessing/PRJEB14954_samples_with_counts_info.txt

# remove unnecessary/tmp files
rm data/1-preprocessing/counts_PRJEB14594_without_header.txt
rm data/1-preprocessing/PRJEB14954_samples_id_with_counts.txt