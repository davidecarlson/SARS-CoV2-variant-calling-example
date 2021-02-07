#!/usr/bin/env python

import argparse
import glob
import gzip
from Bio import SeqIO
import multiprocessing as mp
from itertools import chain
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
		print("SARS-CoV2 assembly size:	{0}".format(length))

# make function to count reads in each fastq

def count_reads(fastq):
	sample = fastq.split('.')[0].split('/')[-1]
	with gzip.open(fastq, "rt") as handle:
			count = 0
			for read in SeqIO.parse(handle, "fastq"):
				count += 1
			return(sample, count)

# make function to get list of length of each read in a fastq file
def read_length(fastq):
	with gzip.open(fastq, "rt") as handle:
		length_list = [len(record.seq) for record in SeqIO.parse(handle, "fastq")]
		return(length_list)

# make function for plotting histogram

def plot_hist(df, x, xlab, png):
	# set style
	sns.set(rc={'figure.figsize':(12,12)})
	sns.set_style("white")
	sns.set_context("paper", font_scale=1.25)
	sns.despine()
	ax1 = sns.displot(df[x], kde=False, bins=8)
	# draw vertical line at mean value
	plt.axvline(stats.mean(df[x]), 0, 0.75, color = "red", linestyle = '--')
	# draw vertical line at median value
	plt.axvline(stats.median(df[x]), 0, 0.75, color = "orange", linestyle = '--')
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
p = mp.Pool(processes=threads)

# use list comprehension to get results
count_results = [p.apply_async(count_reads, [fastq])
        for fastq in fastqs]

p.close()
p.join() # Wait for all child processes to close.
readCount = [p.get() for p in count_results]

#make a dataframe
read_count_df = pd.DataFrame(readCount, columns = ['Sample', 'Read Count'])

#print("The average number of paired end reads per sample is {0} reads.".format(stats.mean(read_count_df['Read Count'])))
print("Average reads per sample:	{0}".format(stats.mean(read_count_df['Read Count'])))
print("Median reads per sample:	{0}".format(stats.median(read_count_df['Read Count'])))
#print("The median number of paired end reads per sample is {0} reads.".format(stats.median(read_count_df['Read Count'])))

# plot the read count histogram

plot_hist(read_count_df, "Read Count", "Read Count", "figs/read_count_hist.png")


# use Multiprocessing to parallelize count of read lengths
p = mp.Pool(processes=threads)


# use list comprehension to get results
length_results = [p.apply_async(read_length, [fastq])
        for fastq in fastqs]

p.close()
p.join() # Wait for all child processes to close.

readLengths = [p.get() for p in length_results]

# flatten the list of lists into a single list
readLength = list(chain(*readLengths))

# make a dataframe
read_length_df = pd.DataFrame(readLength, columns = ['Read Length'])


print("Average read length:	{0}".format(round(stats.mean(read_length_df['Read Length']),2)))
print("Median read length:	{0}".format(stats.median(read_length_df['Read Length'])))

# plot the read length histogram

plot_hist(read_length_df, "Read Length", "Read Length", "figs/read_length_hist.png")