#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=10G,highp
#$ -N indelFiltering
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-323


# Filtering INDELS from GATK:
# start with my variants.vcf : 
# usage qsub script.sh <date you called genotypes in format YYYYMMDD> 
MEM=5G # memory requirement 
# set variables
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load vcftools
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference!
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 
rundate=20171206
PREFIX=PteBra_1
VCF_DIR=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter
vcf=${PREFIX}.raw_variants.${SGE_TASK_ID}.${rundate}.HQsites.Only.1.Variants.vcf.gz

# pull out indels from initial variants file: 

java -Xmx$MEM -jar ${GATK} \
-T SelectVariants \
-R $REFERENCE \
-V ${VCF_DIR}/variants/$vcf \
-selectType INDEL \
-o ${VCF_DIR}/variants/${vcf%.1.Variants.vcf.gz}.i-2.Indels.Variants.vcf.gz

# hard filter the indels using GATK recommendations 
# https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set

java -Xmx$MEM -jar ${GATK} \
-T VariantFiltration \
-R $REFERENCE \
-V ${VCF_DIR}/variants/${vcf%.1.Variants.vcf.gz}.i-2.Indels.Variants.vcf.gz \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filterName "gatk_indel_filter" \
-o ${VCF_DIR}/variants/${vcf%.1.Variants.vcf.gz}.i-3.HF.Indels.Variants.vcf.gz 
    
    
# Want to have PASS only kept:
java -Xmx$MEM -jar ${GATK} \
-R ${REFERENCE} \
-T SelectVariants \
-V ${VCF_DIR}/variants/${vcf%.1.Variants.vcf.gz}.i-3.HF.Indels.Variants.vcf.gz \
-ef \
-o ${VCF_DIR}/variants/${vcf%.1.Variants.vcf.gz}.i-4.PASS-ONLY.HF.Indels.Variants.vcf.gz 

sleep 5m
