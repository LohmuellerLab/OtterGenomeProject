#! /bin/bash
#$ -cwd
#$ -l h_rt=12:00:00,h_data=12G,highp,arch=intel*
#$ -o /u/scratch2/a/ab08028/otters/reports
#$ -e /u/scratch2/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028

# Usage: qsub -N BamFltr01 [script directory here]/run_step7_RemoveBadReads.sh IN.bam

# Optional step to remove bad reads prior to BQSR.
# I use this because my bam files are extremely large, which causes a lot of problems.
# As with MarkDuplicates, I stagger the jobs with a wrapper script due to the IO problems 
# of trying to write many large bam files to disk at once.

source /u/local/Modules/default/init/modules.sh
module load samtools/1.3.1

IN_DIR=/u/scratch2/a/ab08028/otters/bams/IndelRealigned
OUT_DIR=/u/scratch2/a/ab08028/otters/bams/ReadsFiltered
mkdir -p ${OUT_DIR}
samtools view -hb -f 2 -F 256 -q 30 ${IN_DIR}/${1} | samtools view -hb -F 1024 > ${OUT_DIR}/${1%.bam}_Filtered.bam
# http://www.htslib.org/doc/samtools.html
# Explain flags: http://broadinstitute.github.io/picard/explain-flags.html
# hb: include header and output in bam format
# f keep those that match code *numbers are a code not a quantitative thing*
# F discard those that match code
# f 2: Only output alignments with all bits set in INT present in the FLAG field. (2
### 2 == keeping only reads that are "properly aligned according to the aligner"; read mapped in proper pair?
# F 256: Do not output alignments with any bits set in INT present in the FLAG field. 
### 256: is discarding "secondary alignments"
# q: Skip alignments with MAPQ smaller than INT 30
# pass to samtools view again, filter out those with bits present in 1024: get rid of PCR duplicates
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111

PICARD=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Picard_2.8.1/picard.jar

java -jar -Xmx8g -Djava.io.tmpdir=$TMPDIR ${PICARD} BuildBamIndex \
VALIDATION_STRINGENCY=LENIENT TMP_DIR=$TMPDIR \
I=${OUT_DIR}/${1%.bam}_Filtered.bam 

