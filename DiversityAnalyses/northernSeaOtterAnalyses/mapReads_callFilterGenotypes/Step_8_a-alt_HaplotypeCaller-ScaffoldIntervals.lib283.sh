#! /bin/bash
#$ -cwd
#$ -l h_rt=150:00:00,h_data=30G,arch=intel*,highp
#$ -N GTgVCF
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-1
### note! you may be tempted to skip arch=intel* to get the job to start faster
# don't!!! it gives a huge speedup
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
# have to make these before
BAM_DIR=$1 # path to /bams dir 
IN_DIR=$BAM_DIR/ReadsFiltered
OUT_DIR=${BAM_DIR%/bams*}/gvcfs/HaplotypeCaller
PREFIX=$2

mkdir -p $OUT_DIR

intervalFiles=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/intervalFiles

java -jar -Xmx20G $GATK \
-T HaplotypeCaller \
-R $REFERENCE \
-ERC BP_RESOLUTION \
-mbq 20 \
-L $intervalFiles/interval_${SGE_TASK_ID}.bed \
-out_mode EMIT_ALL_SITES \
-I $IN_DIR/${PREFIX}_Aligned_MarkDup_IndelRealigned.JustOneLib.alt.283_Filtered.bam \
-o $OUT_DIR/${PREFIX}.interval_${SGE_TASK_ID}.hapcaller.g.vcf.gz

sleep 5m
