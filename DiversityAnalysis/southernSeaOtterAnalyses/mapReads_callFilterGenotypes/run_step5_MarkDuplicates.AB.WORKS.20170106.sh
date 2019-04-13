#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=32G,highp
#$ -o /u/scratch/a/ab08028/otters/reports
#$ -e /u/scratch/a/ab08028/otters/reports
#$ -m bea
#$ -M ab08028

# Usage: qsub -N MarkDup01 [script directory here]/run_step5_MarkDuplicates.sh IN_lane1.bam IN_lane2.bam OUT.bam

# Note: highly memory-intensive job.
# Running as an array job when bam files are large NOT recommended! Causes IO errors.
# I use a wrapper script to do the submissions one at a time (see submit_run_step5_MarkDuplicates.sh), 
# staggering submissions by a few hours.

# This script as written takes two input bam files to be merged, so edit as needed.

source /u/local/Modules/default/init/modules.sh
module load java

PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar

TEMP_DIR=$TMPDIR
BAM_DIR=/u/scratch/a/ab08028/otters/bams/
mkdir -p $BAM_DIR
BAM1=$1
BAM2=$2
BAM_OUT=$3

cd $BAM_DIR
java -Xmx25G -Djava.io.tmpdir=$TEMP_DIR -jar $PICARD MarkDuplicates \
INPUT=$BAM_DIR/$BAM1 \
INPUT=$BAM_DIR/$BAM2 \
OUTPUT=$BAM_DIR/$BAM_OUT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
METRICS_FILE=$BAM_DIR/${BAM_OUT}_metrics.txt \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
CREATE_INDEX=true \
TMP_DIR=$TEMP_DIR \
MAX_RECORDS_IN_RAM=150000 
