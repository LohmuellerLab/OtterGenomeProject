#! /bin/bash
#$ -cwd
#$ -l h_rt=40:00:00,h_data=10G,highp
#$ -N variantFiltering
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-10
#usage : script <genotype-rundate> <sample size in individuals>
rundate=20171206 # date you called genotypes
ss=1  # POTENTIAL FOR ERRORS HERE; CHECK ss value carefully!!! 

# usage qsub script.sh <date you called genotypes in format YYYYMMDD> 
MEM=5G # memory requirement 
# set variables
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load vcftools
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 

VCF_DIR=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter
PREFIX=PteBra_1
vcf=${PREFIX}.raw_variants.${SGE_TASK_ID}.${rundate}.vcf.gz
# base name but uses HQ Only down below

mkdir -p ${VCF_DIR}/variants

# codes: HF = hard filter; BiSNP = only biallelic SNPs; rmClust = no clustered SNPs (3/10); 

## make this work for my files 
############################################
### Select variants by excluding nonvariants
############################################


java -Xmx$MEM -jar ${GATK} \
    -R ${REFERENCE} \
	-T SelectVariants \
	--excludeNonVariants \
	--removeUnusedAlternates \
	-V ${VCF_DIR}/HQSitesOnly/${vcf%.vcf.gz}.HQsites.Only.rmDotGenotypes.vcf.gz \
	-o ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.1.Variants.vcf.gz
# 20171016: added remove unused alternates to get rid of NON_REF issue if it exists
###############
### Hard filter
###############

java -Xmx$MEM -jar ${GATK} \
    -R ${REFERENCE} \
    -T VariantFiltration \
    -V ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.1.Variants.vcf.gz \
    -filter "QD < 2.0" -filterName "lowQD" \
	-filter "FS > 60.0" -filterName "highFS" \
	-filter "MQ < 40.0" -filterName "lowMQ" \
	-filter "MQRankSum < -12.5" -filterName "lowMQRankSum" \
	-filter "ReadPosRankSum < -8.0" -filterName "lowReadPosRankSum" \
    -o ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.2.HF.Variants.vcf.gz


###################################################
###### Select variants
###### Select biallelic SNP
###### Only retain sites that are annotated as PASS
###################################################

java -Xmx$MEM -jar ${GATK} \
        -R ${REFERENCE} \
        -T SelectVariants \
        -V ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.2.HF.Variants.vcf.gz \
        -selectType SNP \
        --restrictAllelesTo BIALLELIC \
        -o ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.3.BiSNP.HF.Variants.vcf.gz

##################################################
###### Remove clustered SNPs (> 3SNPs within 10bp)
##################################################

java -Xmx$MEM -jar ${GATK} \
        -R ${REFERENCE} \
        -T VariantFiltration \
        -V ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.3.BiSNP.HF.Variants.vcf.gz \
        --clusterWindowSize 10 --clusterSize 3 \
        -o ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.4.rmClust.BiSNP.HF.Variants.vcf.gz


########################
###### Retain PASS sites
########################

java -Xmx$MEM -jar ${GATK} \
-R ${REFERENCE} \
-T SelectVariants \
-V ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.4.rmClust.BiSNP.HF.Variants.vcf.gz \
-ef \
-o ${VCF_DIR}/variants/${vcf%.vcf.gz}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.vcf.gz


sleep 5m
