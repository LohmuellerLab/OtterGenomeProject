#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=500M
#$ -N subFq2sam
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m bea
#$ -M ab08028

# this script will submit run_step2_FastqToSam.AB.sh that converts fastq to uBAM file format

SCRIPTDIR=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/northern_sea_otter_mapping/


QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub

# January 2, 2019:  northern sea otter (enhydra lutris kenyoni) reads from SRA: 
# two lanes of HiSeq X 

# lane 2: https://www.ncbi.nlm.nih.gov/sra/SRX2967283[accn]
# don't use: # lane 1: https://www.ncbi.nlm.nih.gov/sra/SRX2967282[accn] # has a lot of low-quality bases.

# job name , script, input fq 1, input fq 2, output bam name, read group ID (sample_libLane), sample name, library name, flow cell identifier, sequencing center (UCB = Berkeley)

$QSUB -N fq2sam01 $SCRIPTDIR/run_step2_FastqToSam.AB.sh SRX2967283_1.fastq.gz  SRX2967283_2.fastq.gz SRX2967283_FastqToSam.bam Elfin_1 Elfin Lib1 NA Vancouver

