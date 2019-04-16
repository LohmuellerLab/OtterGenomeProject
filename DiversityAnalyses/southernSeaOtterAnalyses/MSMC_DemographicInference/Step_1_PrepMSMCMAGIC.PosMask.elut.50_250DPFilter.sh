#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=10G,highp
#$ -N msmcMagicPrep
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-220
### Already did the masking in Step 11. Can go straight to prep
# for mask, for now just going to use my callable sites. Not even filtering snps yet. (safe? hm. we shall try it a few ways)
# MASKS:
# coverage mask can be my callable sites masks (50-150 mean DP)
# Then snps mask : find complement of sites that are called as variants, but don't end up in my final file, and neg-mask those

# pay attention to the diff between MASK (keeps what's in the mask) and negativemask (removes what's in the mask)
export PATH=$PATH:/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load python/3.4
module load perl
module load vcftools
module load samtools
module load bedtools
PREFIX=01_Elut_CA_Gidget
msmc_tools=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc-tools
magic=/u/home/a/ab08028/klohmueldata/annabel_data/bin/magic
msmc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc
#GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference!
#REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 

# date you called genotypes:
date=20171006
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter

i=${SGE_TASK_ID} # chunk


# THE vcf file has all sites outside of coverage range and indels and genotypes that have a dot (./.) and for bad variants removed

# but you still need to provide a positive mask so msmc knows what sites to consider 
# --mask gives a positive mask (sites to KEEP)
# --negative_mask gives a neg. mask (sites to discard)
# PREP FOR MSMC: 
mkdir -p $wd/msmcAnalysis
mkdir -p $wd/msmcAnalysis/inputFiles
posMask=$wd/FinalGoodSitesPosMasks_forMSMC/chunk_${i}_PosMask.merged.bed.gz
inputvcf=$wd/HQSitesOnly-MaskedBadVariantsIndels/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz
$msmc_tools/generate_multihetsep.py $inputvcf --mask $posMask > $wd/msmcAnalysis/inputFiles/chunk_${i}_postMultiHetSep.txt

# optional: PREP FOR MAGIC:
#mkdir -p $wd/magicAnalysis
#mkdir -p $wd/magicAnalysis/magicCoverageFiles
#mkdir -p $wd/magicAnalysis/magicPrep
# combo prep:
## THIS isn't working right; mask file may need to be unzipped? skipping for now. 
# $magic/combo_prep.py $inputvcf --masks $posMask --coverfile $wd/magicAnalysis/magicCoverageFiles/chunk_${i}_cover.dist.txt $wd/magicAnalysis/magicPrep/chunk_${i}_postComboPrep.txt > $wd/magicAnalysis/magicPrep/chunk_${i}_postComboPrep.txt
# windower (with default tbl setting)
# $magic/windower.py --coverfile $wd/magicAnalysis/magicCoverageFiles/chunk_${i}_cover.dist.txt $wd/magicAnalysis/magicPrep/chunk_${i}_postComboPrep.txt 
# this outputs in the same directory, with "counts" as suffix. 

sleep 5m
