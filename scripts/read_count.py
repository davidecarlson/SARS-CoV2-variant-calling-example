#!/usr/bin/env python

import argparse
import glob
import gzip
from Bio import SeqIO
import multiprocessing as mp
import pandas as pd
import seaborn as sns
import statistics  as stats
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Get some basic stats about the SARS-CoV2 data set")

parser.add_argument('--threads', type=int, required=False, default=1, help='Number of threads to use. Default is 1', action='store')

args=parser.parse_args()

threads = args.threads

# make function to calculat the size of the reference genome

def genome_size(fasta):
	for record in SeqIO.parse(fasta, "fasta"):
		length = len(record.seq)
		#print(record.id + "\t" + str(len(record.seq)))
		print("The size of the SARS-CoV2 reference assembly is {0} base pairs".format(length))

# make function to count reads in each fastq

def count_reads(fastq):
	sample = fastq.split('.')[0].split('/')[-1]
	with gzip.open(fastq, "rt") as handle:
			count = 0
			for read in SeqIO.parse(handle, "fastq"):
				count += 1
			return(sample, count)

# make function for plotting histogram

def plot_hist(df, xlab, png):
	# set style
	sns.set(rc={'figure.figsize':(12,12)})
	sns.set_style("white")
	sns.set_context("paper", font_scale=1.25)
	sns.despine()
	ax1 = sns.displot(df['Read Count'], kde=False)
	plt.axvline(stats.mean(df['Read Count']), 0, 0.75, color = "red", linestyle = '--')
	plt.xticks(rotation=45)
	plt.xlabel(xlab)
	plt.savefig(png, bbox_inches='tight')
	plt.clf()

# get the genome assembly size

genome_size('assembly/GCF_009858895.2_ASM985889v3_genomic.fna')

# specify directory with fastq files
readDir = 'fastqs/'

fastqs = glob.glob(readDir + '*1.fastq.gz')


# use Multiprocessing to parallelize counting reads
p = mp.Pool(processes=25)

# use list comprehension to get results
results = [p.apply_async(count_reads, [fastq])
        for fastq in fastqs]

p.close()
p.join() # Wait for all child processes to close.
output = [p.get() for p in results]

read_df = pd.DataFrame(output, columns = ['Sample', 'Read Count'])

# plot the read count histogram

plot_hist(read_df, "Read Count", "figs/read_count_hist.png")

