#! /bin/bash
#$ -wd /u/flashscratch/a/ab08028/otters/bams
#$ -l h_rt=56:00:00,h_data=10G,highp
#$ -pe shared 3
#$ -m bea
#$ -M ab08028
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports

source /u/local/Modules/default/init/modules.sh
module load java
module load bwa/0.7.7
PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar


REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

BAM_DIR=/u/flashscratch/a/ab08028/otters/bams

FILENAME=${1}

TEMP_DIR=/u/flashscratch/a/ab08028/temp-genome
mkdir -p $TEMP_DIR

cd $BAM_DIR

set -o pipefail

java -Xmx8G -jar -Djava.io.tmpdir=$TEMP_DIR \
$PICARD SamToFastq \
I=$BAM_DIR/${FILENAME}_MarkIlluminaAdapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=$TEMP_DIR 2>>./"Process_"${FILENAME}"_SamToFastq.txt" | \
bwa mem -M -t 3 -p $REFERENCE /dev/stdin 2>>./"Process_"${FILENAME}"_BwaMem.txt" | \
java -Xmx8G -jar -Djava.io.tmpdir=$TEMP_DIR \
$PICARD MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=$BAM_DIR/${FILENAME}_FastqToSam.bam \
OUTPUT=$BAM_DIR/${FILENAME}_Aligned.bam \
R=$REFERENCE CREATE_INDEX=true \
ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=$TEMP_DIR 2>>./"Process_"${FILENAME}"_MergeBam.txt"
