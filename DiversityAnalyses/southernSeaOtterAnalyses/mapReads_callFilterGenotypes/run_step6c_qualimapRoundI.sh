#! /bin/bash
#$ -cwd
#$ -l h_rt=60:00:00,h_data=2G,highp
#$ -pe shared 15
#$ -N qualimap
#$ -o output.txt
#$ -e error.txt
#$ -m abe
#$ -M ab08028

qmap=/u/home/a/ab08028/klohmueldata/annabel_data/bin/qualimap_v2.2.1/qualimap

source /u/local/Modules/default/init/modules.sh
module load R
module load java
INDIR=/u/scratch2/a/ab08028/otters/bams/IndelRealigned

$qmap bamqc --java-mem-size=24G -nt 15 -bam $INDIR/01_Elut_CA_Gidget_Aligned_MarkDup_IndelRealigned.bam --outfile qualimap.out.pdf

# outfile will end up in directory with input bam
