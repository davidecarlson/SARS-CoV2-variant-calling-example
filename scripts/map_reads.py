#!/usr/bin/env python

import glob
import os
import subprocess
import multiprocessing as mp
import argparse
from itertools import repeat

parser = argparse.ArgumentParser(description="Get some basic stats about the SARS-CoV2 data set")

parser.add_argument('--threads', type=int, required=False, default=1, help='Number of threads to use. Default is 1', action='store')

args=parser.parse_args()

threads = args.threads

# run bwa mem, samtools, and picard MarkDup on a series of fastq files all found in the same folder
# assumes a fastq name structure of <sample id>_< 1 or 2>.fastq.gz

# path to reference genome
reference = 'assembly/GCF_009858895.2_ASM985889v3_genomic.fna' # must be indexed already

# path to directory with reads

read_dir = 'fastqs'

fwd_reads = glob.glob(read_dir + '/*_1.fastq.gz')

sample_list = [read.split('_')[0].split('/')[-1] for read in fwd_reads]

# create output directories if they don't exist

bams = 'bams'
logs = 'logs'
stats = 'stats'

outdirs = [bams, logs, stats]

for dir in outdirs:
	if not os.path.exists(dir):
		os.makedirs(dir)

def map_reads(sample,in_dir, bam_dir, log_dir, stats_dir, ref, threads):
	# create variables for running bwa
	fwd = in_dir + '/' + sample + '_1.fastq.gz'
	rev = in_dir + '/' + sample + '_2.fastq.gz'

	# create names for output files:
	bamfile = bams + '/' + sample + '.bam'
	sortedbam = bams + '/' + sample + '.sorted.bam'
	markdupbam = bams + '/' + sample + '.sorted.markdup.bam'
	metrics = stats + '/' + sample + '.metrics.txt '
	index = bams + '/' + sample + '.sorted.markdup.bam.bai'

	# check if bam file exists already and looks good
	bamcheck = subprocess.run(['samtools', 'quickcheck', markdupbam], capture_output=True, text=True)

	return_code = bamcheck.returncode
	if return_code == 0:
		print("Read mapping already completed for sample {0}!".format(sample))
	elif return_code > 0:

		# Map the reads to the reference
		# add a read group
		rg=r"'@RG\tID:"+ sample+r"\tSM:" + sample +r"\tLB:lib1\tPL:ILLUMINA'"
		bwa_command = 'bwa mem -R ' + rg + ' -t ' + str(threads) + ' ' + ref + ' ' + fwd + ' ' + rev + ' | samtools view -bS  -o ' + bamfile
		# create logfile
		bwalogfile = open(logs + '/' + sample + ".bwa.log", "w")
		# run bwa to map reads and get sam file
		print("Mapping reads to reference for sample {0}!".format(sample))
		subprocess.call(bwa_command, stderr=bwalogfile, stdout=bwalogfile, shell=True)

		# sort the bam file by position
		print("Sorting bam for sample {0}!".format(sample))
		samtools_command = 'samtools sort --threads ' + str(threads) + ' -o ' + sortedbam + ' ' + bamfile
		subprocess.call(samtools_command, stderr=subprocess.STDOUT, shell=True)
		# remove unsorted bam
		os.remove(bamfile)

		# mark pcr and optical duplicates
		print("Marking duplicates for sample {0}!".format(sample))
		markdup_command = 'picard MarkDuplicates --INPUT ' + sortedbam + ' --METRICS_FILE ' + metrics + ' --OUTPUT ' + markdupbam
		# create logfile
		duplogfile = open(logs + '/' + sample + ".markdup.log", "w")
		subprocess.call(markdup_command, stderr=duplogfile, stdout=duplogfile, shell=True)
		os.remove(sortedbam)

		# index the final bam files
		index_command = 'samtools index ' + markdupbam + ' ' + index
		subprocess.call(index_command, stderr=subprocess.STDOUT, shell=True)


# only using 1 thread per sample - this isn't efficient, but it's just an example and we're running (optionally) more than one sample a time

with mp.Pool(processes=threads) as pool:
	results = pool.starmap_async(map_reads, zip(sample_list, repeat(read_dir), repeat(bams), repeat(logs), repeat(stats), repeat(reference), repeat(1)))
	results.get()

pool.close()
pool.join() # Wait for all child processes to close.

print("Finished mapping reads for all samples!")