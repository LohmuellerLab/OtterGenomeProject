#! /bin/bash
#$ -cwd
#$ -l h_rt=12:00:00,h_data=21G,arch=intel*
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m bea
#$ -M ab08028
#$ -N fastq2Sam

# Usage: qsub -N fq2sam01 IN_Read1.fastq.gz IN_Read2.fastq.gz OUT_FastqToSam.bam [RGID] [RGSM] [RGLB] [RGPL] [RGCN]
## This script will convert fastq file to uBAM file. You submit it to Hoffman using submit_run_step2


# RGID is unique identifier (you come up with, e.g. gidget_1a for gidget, library 1, lane a)
# RGSM is sample name (e.g. GIDGET)
# RGLB is library name (e.g. Lib1)
# RGPL is platform: illumina
# RGCN is sequencing center: UCB
# See the wrapper script, submit_run_FastqToSam.sh, for job submission example.

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111 # need to be java 1.8


PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar

READ_DIR=/u/flashscratch/a/ab08028/otters/fastqs
BAM_OUTDIR=/u/flashscratch/a/ab08028/otters/bams
TEMP_DIR=/u/flashscratch/a/ab08028/temp

java -Xmx16G -jar -Djava.io.tmpdir=$TEMP_DIR \
$PICARD FastqToSam \
FASTQ=$READ_DIR/$1 \
FASTQ2=$READ_DIR/$2 \
OUTPUT=$BAM_OUTDIR/$3 \
READ_GROUP_NAME=$4 \
SAMPLE_NAME=$5 \
LIBRARY_NAME=$6 \
PLATFORM_UNIT=$7 \
PLATFORM=illumina \
SEQUENCING_CENTER=$8 \
TMP_DIR=$TEMP_DIR

