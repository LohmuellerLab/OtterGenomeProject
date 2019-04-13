#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=15G,highp
#$ -N maskHQFiles
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-323

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111


GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

# ferret reference!
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 

# want to remove the failed variants (negmask1? 3?) and the yucky sites from HQ Sites vcf
rundate=20171006
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter
mkdir -p $wd/HQSitesOnly-MaskedBadVariantsIndels
PREFIX=01_Elut_CA_Gidget
vcf=${PREFIX}.raw_variants.${SGE_TASK_ID}.${rundate}.vcf.gz
i=${SGE_TASK_ID}
# Exclude sites that failed variant calling (Step 11)
# and Exclude any sites where the reference allele is multinucleotide 
# See Step 11 b for how masks were made
# negmask2 are variants with ./. genotype; were already masked in Step 10e-ii
# negmask1 : these are variants that didn't end up in final file (includes indels)
negmask1=$wd/variants/badVarsNegativeMasksForMSMC/chunk_${i}_badVars.0based.bed
# negmask3: these are weird variants (deletions-ish) where ref allele is multinucleotide 
negmask3=$wd/multiNucRefAlleleBedCoords/chunk_${i}.multinuc.RefAllele.0based.bed
# XL excludes sites:
java -jar -Xmx10G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-XL ${negmask1} \
-XL ${negmask3} \
-V $wd/HQSitesOnly/${vcf%.vcf.gz}.HQsites.Only.rmDotGenotypes.vcf.gz   \
-o $wd/HQSitesOnly-MaskedBadVariantsIndels/${vcf%.vcf.gz}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz

sleep 5m
