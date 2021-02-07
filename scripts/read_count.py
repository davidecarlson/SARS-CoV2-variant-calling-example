#!/usr/bin/env python

import glob
import gzip
from Bio import SeqIO
import multiprocessing as mp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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
	sns.set_context("paper", font_scale=1.75)
	sns.despine()
	ax1 = sns.displot(df['Read Count'], kde=False)
	plt.axvline(mean(df['Read Count'], 0, 0.75))
	plt.xlabel(xlab)
	plt.savefig(png, bbox_inches='tight')

# specify directory with fastq files
readDir = '~/coronavirus_example/fastqs/'

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

plot_hist(read_df, "Read Count", "~/coronavirus_example/figs/read_count_hist.png")

