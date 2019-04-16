# script to count number of genotypes called and number of heterozygotes per sample
# input file is a VCF file that has been filtered
# supply input file, window size, step size (don't need chr number of iterator ,will get it from vcf file)
# what happens if vcf has more than 1 chr? we'll see.
# example: python ./SlidingWindowHet.py input.vcf.gz 100000 10000 fai_file (1-> total # of chrom?)

import sys
import pysam
import os
import gzip
import csv
import argparse

def parse_args():
	"""
	Parse command-line arguments
	"""
	parser = argparse.ArgumentParser(description="This script computes sliding window heterozygosity.")

	parser.add_argument(
            "--vcf", required=True,
            help="REQUIRED. input vcf file, this must have passed through Step 11 so *NO* bad sites remain in the file")
	parser.add_argument(
			"--window_size", required=True,
			help="REQUIRED. window size ")
	parser.add_argument("--step_size", required=True,
			help="REQUIRED. Name of output file.")
	parser.add_argument("--fai", required=True,
			help="REQUIRED. fai file of genome")
	args = parser.parse_args()
	return args
args=parse_args()
# ab: I have changed this so you use a VCF that has had all not-pass sites removed. HQ callable sites only, bad variants and indels removed
# open input file (gzipped VCF file), make sure the VCF file is indexed (if not, create index)
filename=args.vcf
VCF = gzip.open(filename, 'r')

if not os.path.exists("%s.tbi" % filename):
	pysam.tabix_index(filename, preset="vcf")
parsevcf = pysam.Tabixfile(filename)


# set variables
window_size = int(args.window_size)
step_size = int(args.step_size)

# get scaffold sizes

scaffs=[]
lengths=[]
with open(args.fai,"r") as faiFile:
    for line in faiFile:
        line = line.split("\t")
        scaffs.append(line[0])
        lengths.append(line[1])
        chromo_size=dict(zip(scaffs,lengths))

# get list of samples from the line that starts with a single #
samples=[]
for line in VCF:
	if line.startswith('##'):
		pass
	else:
		for i in line.split()[9:]: samples.append(i)
		break


# get first and last positions in chromosome
# this looks at the first non # line:
for line in VCF:
	if line[0] != '#':
		chromo=str(line.strip().split()[0])
		start_pos = int(line.strip().split()[1])
		end_pos = int(chromo_size[chromo]) # get chr name from the first field of the vcf , look up its length in chromo_size
		break


# create output file
output = open(filename + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
output.write('chromo\twindow_start\twindow_end\tchrLength\tsites_total\tsites_unmasked\tsites_passing\tcalls_%s\thets_%s\n' % ('\tcalls_'.join(samples), '\thets_'.join(samples)) )


# Fetch a region, ignore sites that fail filters, tally total calls and heterozygotes
def snp_cal(chromo,window_start,window_end):

	rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chromo, window_start, window_end), parser=pysam.asTuple()))

	sites_total,sites_unmasked,sites_passing=0,0,0
	calls=[0]*len(samples)
	hets=[0]*len(samples)

	for line in rows:
		sites_total+=1
		#if "CpGRep" in line[6]: continue
		sites_unmasked+=1
		#if "FAIL" in line[6]: continue
		#if "WARN" in line[6]: continue
		sites_passing+=1
		for i in range(0,len(samples)):
			GT=line[i+9]
			if GT[:1]!='.': calls[i]+=1
			if GT[:3]=='0/1': hets[i]+=1

	output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chromo,window_start,window_end,chromo_size[chromo],sites_total,sites_unmasked,sites_passing,'\t'.join(map(str,calls)),'\t'.join(map(str,hets))) )


# initialize window start and end coordinates
window_start = start_pos
window_end = start_pos+window_size-1


# calculate stats for window, update window start and end positions, repeat to end of chromosome
while window_end <= end_pos:

	if window_end < end_pos:
		snp_cal(chromo,window_start,window_end)

		window_start = window_start + step_size
		window_end = window_start + window_size - 1

	else:
		snp_cal(chromo,window_start,window_end)
		break

else:
	window_end = end_pos
	snp_cal(chromo,window_start,window_end)


# close files and exit
VCF.close()
output.close()

exit()
