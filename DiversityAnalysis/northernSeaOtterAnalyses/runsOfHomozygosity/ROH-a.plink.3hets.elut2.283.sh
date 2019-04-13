#! /bin/bash
#$ -cwd
#$ -l h_rt=02:00:00,h_data=10G,highp
#$ -N plinkROH
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-220

 
# vcftools LROH requires >1 individual
# plink doesn't

source /u/local/Modules/default/init/modules.sh
module load vcftools
module load plink

# need to get chr name from file
i=$SGE_TASK_ID
date=20190215_nso_lib283only # date of vcf calling
PREFIX=02_Elut_SEAK_Elfin
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/genome-stats/ROH
mkdir -p $wd
vcf=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/HQSitesOnly-MaskedBadVariantsIndels/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz

chr=`zcat $vcf | grep -v "#" | head -n 1 | awk '{print $1}'`


# convert to ped/map format. note that you lose chromosome info - that's a pain. 
plinkindir=$wd/plinkInputFiles
plinkoutdir=$wd/plinkOutputFiles/
mkdir -p $plinkoutdir
mkdir -p $plinkindir
vcftools --gzvcf $vcf --plink --chr $chr --out $plinkindir/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.Plink

# use plink:
# defaults (also cite MacQuillan 2008)
# sliding window is 5000kb, 50 snps
# window allows 5 missing calls, 3 het call
# 1000kb min length and >100 variants for a ROH --> CHANGING this to 500kb and >50 variants ala MacQuillan 
# at least 1 variant per 50kb
# http://zzz.bwh.harvard.edu/plink/ibdibs.shtml
# http://www.sciencedirect.com/science/article/pii/S000292970800445X


plink \
--file $plinkindir/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.Plink \
--homozyg \
--homozyg-window-het 3 \
--homozyg-window-missing 5 \
--homozyg-kb 500 \
--homozyg-snp 50 \
--out $plinkoutdir/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.Plink.out


sleep 5m
