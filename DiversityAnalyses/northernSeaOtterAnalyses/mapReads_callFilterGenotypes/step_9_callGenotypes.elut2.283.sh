#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=32G,highp
#$ -N jointGeno
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-323


# note: sometimes finishes faster than anticipated (a few hours) and seems like it hasn't finished, but if file isn't truncated
# it is apparently okay; schedule is just unreliable https://gatkforums.broadinstitute.org/gatk/discussion/5329/joint-genotyping


rundate=20190215_nso_lib283only
wd=/u/flashscratch/a/ab08028/otters/
PREFIX=02_Elut_SEAK_Elfin
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load samtools
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta


OUT_DIR=$wd/vcfs/vcf_${rundate}/raw_variants
gVCF_DIR=/u/flashscratch/a/ab08028/otters/gvcfs/HaplotypeCaller

mkdir -p ${OUT_DIR}


java -Xmx15G -jar $GATK \
  -T GenotypeGVCFs \
  -R $REFERENCE \
  -V ${gVCF_DIR}/${PREFIX}.interval_${SGE_TASK_ID}.hapcaller.g.vcf.gz \
  --includeNonVariantSites \
  -o ${OUT_DIR}/${PREFIX}.raw_variants.${SGE_TASK_ID}.${rundate}.vcf.gz
  
# check for completion!
sleep 10m
