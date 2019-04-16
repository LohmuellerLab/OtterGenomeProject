#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=28G,arch=intel*,highp
#$ -N getCallableSites
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-323

## note: the output of this will be very large *lots of warn messages* -- gzip the reports afterward

##########*** NOTE manually put vcf files into a directory called raw_variants * #############

source /u/local/Modules/default/init/modules.sh
module load python/2.7
module load java
module load vcftools
module load bedtools

#usage : script <genotype-rundate> <sample size in individuals> 


rundate=20190215_nso_lib283only
 # date you called genotypes if doing this on a different day as 
ss=1  # POTENTIAL FOR ERRORS HERE; CHECK ss value carefully!!! 
PREFIX=02_Elut_SEAK_Elfin

## get this from steps 10a-d (using DP dist)
# elut2 stats:
# mean = 39
# 50% min
minDP=19 # 
# 250% max
maxDP=97 # 

GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

# ferret reference
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

SCRIPT_DIR=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/northern_sea_otter_mapping


IN_DIR=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}/raw_variants # initial in dir is the where the raw gts are
VCF_DIR=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250_Filter
mkdir -p ${VCF_DIR}
# for now: have outdir be sandbox
vcf=${PREFIX}.raw_variants.${SGE_TASK_ID}.${rundate}.vcf.gz  

mkdir -p ${VCF_DIR}/HQSitesMarked
mkdir -p ${VCF_DIR}/HQSitesBedCoords
mkdir -p ${VCF_DIR}/HQSitesOnly-stillbadgenotypes
# get HQ sites: 
java -jar -Xmx16G ${GATK} \
-T VariantFiltration \
-R ${REFERENCE} \
-V ${IN_DIR}/$vcf \
--filterExpression "AN < $((2*$ss))" --filterName "MissingCalledGenotypes" \
--filterExpression "DP < ${minDP} || DP > ${maxDP} " --filterName "DP_Filter_20x_250x_new" \
-o ${VCF_DIR}/HQSitesMarked/${vcf%.vcf.gz}.HQsites.Marked.vcf.gz

# get coordinates of HQ sites using Tanya's script:

obtainHQ=${SCRIPT_DIR}/obtain_high_quality_coordinates.py # Tanya's script  # scriptdir ok
python $obtainHQ --VCF ${VCF_DIR}/HQSitesMarked/${vcf%.vcf.gz}.HQsites.Marked.vcf.gz --outfile ${VCF_DIR}/HQSitesBedCoords/chunk_${SGE_TASK_ID}_HQsites.coords.0based.bed
# you can then merge the bed file to make it more efficient (works either way)

# Filter the VCF to only have those coodinates; make sure it's sorted first
# note that you have to sort a bed file before it can be merged! Tanya's script has it pre-sorted
# so you can go sraight to merge (way faster)
bedtools merge -i ${VCF_DIR}/HQSitesBedCoords/chunk_${SGE_TASK_ID}_HQsites.coords.0based.bed  > ${VCF_DIR}/HQSitesBedCoords/chunk_${SGE_TASK_ID}_HQsites.coords.0based.merged.bed

Bed=${VCF_DIR}/HQSitesBedCoords/chunk_${SGE_TASK_ID}_HQsites.coords.0based.merged.bed

# pull out only the HQ sites:
java -jar -Xmx16G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-L ${Bed} \
-V ${VCF_DIR}/HQSitesMarked/${vcf%.vcf.gz}.HQsites.Marked.vcf.gz   \
-o ${VCF_DIR}/HQSitesOnly-stillbadgenotypes/${vcf%.vcf.gz}.HQsites.Only.vcf.gz


sleep 5m
