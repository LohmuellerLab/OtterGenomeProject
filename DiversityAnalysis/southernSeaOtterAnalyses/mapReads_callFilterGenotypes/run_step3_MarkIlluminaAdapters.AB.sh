#! /bin/bash
#$ -wd /u/scratch/a/ab08028/otters/
#$ -l h_rt=12:00:00,h_data=21G,arch=intel*,highp
#$ -o /u/scratch/a/ab08028/otters/reports
#$ -e /u/scratch/a/ab08028/otters/reports
#$ -m bea
#$ -M ab08028

# Usage: qsub -N jobname run_step3_MarkIlluminaAdapters.AB.sh INPUT_FastqToSam.bam


source /u/local/Modules/default/init/modules.sh
module load java

cd /u/scratch/a/ab08028/otters/bams

FILENAME=$1
PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar

TEMP_DIR=/u/scratch/a/ab08028/temp

java -Xmx16G -jar -Djava.io.tmpdir=$TEMP_DIR \
$PICARD \
MarkIlluminaAdapters \
I=$FILENAME \
O=${FILENAME%_FastqToSam.bam}_MarkIlluminaAdapters.bam \
M=${FILENAME%_FastqToSam.bam}_MarkIlluminaAdapters.bam_metrics.txt \
TMP_DIR=$TEMP_DIR
