# SARS-CoV2-variant-calling-example

This repository includes code and instructions to perform an example variant calling analysis on SARS-CoV2 data. All data used in the analysis are publicly available.

####Step 1
Download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

####Step 2
Clone this repository:
```
git clone https://github.com/davidecarlson/SARS-CoV2-variant-calling-example.git
```

####Step 3
Create a new conda environment with the software needed in this analysis using the provided environment file:
```
conda env create -f environment.yaml
````

####Step 4
Download the SARS-CoV2 reference genome (Accession NC_045512.2) in fasta format from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)