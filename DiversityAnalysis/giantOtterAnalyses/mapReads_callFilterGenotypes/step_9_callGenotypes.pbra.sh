#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=20G
#$ -N jointGeno
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-323

rundate=20171206
wd=/u/flashscratch/a/ab08028/otters/
PREFIX=PteBra_1
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load samtools
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 


OUT_DIR=$wd/vcfs/vcf_${rundate}
gVCF_DIR=/u/flashscratch/a/ab08028/otters/gvcfs/HaplotypeCaller

mkdir -p ${OUT_DIR}
bgzip=~/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

# only genotyping 1 vcf 

java -Xmx15G -jar $GATK \
  -T GenotypeGVCFs \
  -R $REFERENCE \
  -V $gVCF_DIR/${PREFIX}.interval_${SGE_TASK_ID}.hapcaller.g.vcf.gz \
  --includeNonVariantSites \
  -o ${OUT_DIR}/${PREFIX}.raw_variants.${SGE_TASK_ID}.${rundate}.vcf.gz
  

sleep 5m
