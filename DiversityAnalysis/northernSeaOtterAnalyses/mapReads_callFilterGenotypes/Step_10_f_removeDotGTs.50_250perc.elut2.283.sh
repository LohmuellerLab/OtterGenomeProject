#! /bin/bash
#$ -cwd
#$ -l h_rt=1:00:00,h_data=5G
#$ -N removeDotGenotypes
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-323

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load bedtools
module load samtools/1.3.1
bgzip=~/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix


# date you called genotypes:
date=20190215_nso_lib283only
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter

i=${SGE_TASK_ID} # chunk
mkdir -p $wd/HQSitesOnly

PREFIX=02_Elut_SEAK_Elfin
zcat $wd/HQSitesOnly-stillbadgenotypes/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.vcf.gz | \
grep -v "\./\." | $bgzip > $wd/HQSitesOnly/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.vcf.gz

$tabix -p vcf $wd/HQSitesOnly/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.vcf.gz

sleep 5m

