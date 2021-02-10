# SARS-CoV2-variant-calling-example

This repository includes code and instructions to perform an example variant calling analysis on SARS-CoV2 data. All data used in the analysis are publicly available.

#### Step 1
Download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

#### Step 2
Clone this repository:
```
git clone https://github.com/davidecarlson/SARS-CoV2-variant-calling-example.git
```

#### Step 3
Create a new conda environment with the software needed in this analysis using the provided environment file:
```
conda env create -f environment.yaml
````

#### Step 4
Activate the environment:
```conda activate covid-variant-example
```

#### Step 5
Download the SARS-CoV2 reference genome (Accession NC_045512.2) in fasta format from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)

#### Step 6
Execute scripts/download_files.sh to download 100 sets of paired-end Illumina reads from BioProject PRJNA656534
```
# read the accesion numbers listed in fastqs/samples.txt and use GNU parallel and fastq-dump to download paired end reads for each sample
# change the value of the argument provided to "--jobs" to some number equal to or less than the number of cores on your machine
./scripts/download_files.sh
```

#### Step 7
Execute scripts/data_stats.py to get some basic stats/plots about the reference genome and fastq files:
```
python scripts/data_stats.py --threads <number of desired threads>
```
If all goes well, the following lines should be printed to the screen:
```
SARS-CoV2 assembly size:        29903
Average reads per sample:       271432.66
Median reads per sample:        176628.0
Average read length:    117.49
Median read length:     136.0
```

#### Step 8

Execute scripts/map_reads.py to map the paired-end reads for each sample against the reference genome:

```
python scripts/map_reads.py --threads <number of samples to map in parallel>
```