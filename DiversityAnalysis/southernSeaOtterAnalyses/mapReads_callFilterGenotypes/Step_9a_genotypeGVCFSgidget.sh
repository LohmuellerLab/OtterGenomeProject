#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=10G,highp,arch=intel*
#$ -N jointGeno
#$ -o /u/scratch2/a/ab08028/otters/reports
#$ -e /u/scratch2/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-323
# this runs very quickly! don't need to give it 5 hrs

rundate=`date +%Y%m%d` 
wd=/u/scratch2/a/ab08028/otters/
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load samtools
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 


OUT_DIR=$wd/vcfs/vcf_${rundate}
gVCF_DIR=/u/scratch2/a/ab08028/otters/gvcfs/HaplotypeCaller

mkdir -p ${OUT_DIR}
bgzip=~/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip


java -jar $GATK \
  -T GenotypeGVCFs \
  -R $REFERENCE \
  -V $gVCF_DIR/01_Elut_CA_Gidget.interval_${SGE_TASK_ID}.hapcaller.g.vcf.gz \
  --includeNonVariantSites \
  -o ${OUT_DIR}/01_Elut_CA_Gidget.raw_variants.${SGE_TASK_ID}.${rundate}.vcf.gz
  

sleep 10m
