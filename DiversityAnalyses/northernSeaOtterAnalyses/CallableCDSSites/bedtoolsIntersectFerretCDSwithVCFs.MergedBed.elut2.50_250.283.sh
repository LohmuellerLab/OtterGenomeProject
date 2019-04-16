#!/bin/bash
#$ -N intersectVCFWithCodingSeq
#$ -l h_data=2G,h_rt=10:00:00
#$ -m bea
#$ -M ab08028
#$ -cwd
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports

### Bedtools intersect
# Get coding sequences from vcf file
source /u/local/Modules/default/init/modules.sh
module load bedtools
date=20190215_nso_lib283only # northern sea otter 
header=02_Elut_SEAK_Elfin
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter

# now going to make gff into a bed file
# RECALL that gff files are 1-based
# but bed files are 0-based non inclusive, so start site needs 1 subtracted, but end site is not inclusive so should stay the same (subtracting one is already part of it)
## only need to do this once! and then can use it for pbra. 
# already did these steps:
#awk '{OFS="\t"; print $1,$4-1,$5}' $gff > ${gff%.gff}.0based.bed
#bed=${gff%.gff}.0based.bed
# needs to be sorted by chrom and then by start 
#bedtools sort -i $bed > ${bed%.bed}.sorted.bed
#bedtools merge -i ${bed%.bed}.sorted.bed > ${bed%.bed}.sorted.merged.bed

bed=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/MusPutFuro1.0.91.cdsOnly.0based.sorted.merged.bed

mkdir -p $wd/cds_vcfs

for i in {1..323}
do
vcf=${header}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz
bedtools intersect -header -a $wd/HQSitesOnly-MaskedBadVariantsIndels/$vcf -b $bed -wa > $wd/cds_vcfs/${header}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf
done
sleep 5m
