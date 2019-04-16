#!/bin/bash
#$ -N intersectVCFWithCodingSeq
#$ -l h_data=2G,h_rt=10:00:00,highp
#$ -m bea
#$ -M ab08028
#$ -cwd
### note for indels you're intersecting the variant sites, whereas SNPs we did ALL sites (won't make a difference,
# still just want to annotate the variant sites anyway)
### Bedtools intersect
# Get coding sequences from vcf file
source /u/local/Modules/default/init/modules.sh
module load bedtools
date=20171006 # is sea otters; 20171206 is giant otter
PREFIX=01_Elut_CA_Gidget
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter
# make a bed file from the ferret gff file of cds sequences
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

mkdir -p $wd/cds_vcfs/indels

for i in {1..323}
do
vcf=${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.vcf.gz
bedtools intersect -header -a $wd/variants/$vcf -b $bed -wa | sort | uniq > $wd/cds_vcfs/indels/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.CDSOnly.fromMergedBed.vcf
## adding in sort | uniq so that sites are only reported once
done
sleep 5m
