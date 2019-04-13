#!/bin/bash
#$ -l h_rt=1:00:00,h_data=2G,highp
#$ -N makePosMaskForMSMC
#$ -cwd
#$ -m bea
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -M ab08028
#$ -t 1-220
# Make a bed file of callable sites and merge it and bgzip it
# first 220 scaffolds
i=${SGE_TASK_ID}
PREFIX=01_Elut_CA_Gidget
rundate=20171006
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter
outdir=$wd/FinalGoodSitesPosMasks_forMSMC
mkdir -p $outdir
source /u/local/Modules/default/init/modules.sh
module load bedtools
# want to make a bed file of all the sites
### use the vcf where you've removed all bad things (./., indels, bad variants, non callable sites, etc.)
# just want a mask with these coordinates 
vcf=HQSitesOnly-MaskedBadVariantsIndels/$PREFIX.raw_variants.${i}.${rundate}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz
zcat $wd/$vcf | grep -v "#" | awk '{OFS="\t";print $1,$2-1,$2}'> $outdir/chunk_${i}_PosMask.bed
bedtools merge -i $outdir/chunk_${i}_PosMask.bed > $outdir/chunk_${i}_PosMask.merged.bed
gzip $outdir/chunk_${i}_PosMask.bed
gzip $outdir/chunk_${i}_PosMask.merged.bed

sleep 5m