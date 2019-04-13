#! /bin/bash
#$ -wd /u/scratch/a/ab08028/otters/bams
#$ -l h_rt=24:00:00,h_data=1G
#$ -N subBamFltr
#$ -o /u/scratch/a/ab08028/otters/reports
#$ -e /u/scratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
# do this after indel realignment (reduces bam size a lot)
# do qualimap before and after

qsub /u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts/run_step7_filterBadReads.optional.sh 01_Elut_CA_Gidget_Aligned_MarkDup_IndelRealigned.bam
