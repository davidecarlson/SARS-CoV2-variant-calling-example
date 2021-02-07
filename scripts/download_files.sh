#!/usr/bin/env bash

# download paired-end fastq files for every sample in BioProject PRJNA656534

FASTQDIR=~/coronavirus_example/fastqs/
SAMPLES=~/coronavirus_example/fastqs/samples.txt

cat $SAMPLES | parallel --jobs 25 --verbose "fastq-dump --accession {} --split-files --gzip --outdir $FASTQDIR"
