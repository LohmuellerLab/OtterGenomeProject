#! /bin/bash
#$ -cwd
#$ -l h_rt=05:00:00,h_data=5G,highp
#$ -N slidingWindowHet
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-220
source /u/local/Modules/default/init/modules.sh
module load python/2.7
# run sliding window het for chunk 1:
date=20171006
PREFIX=01_Elut_CA_Gidget
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter
i=${SGE_TASK_ID}
window=1000000
step=20000
vcf=$wd/HQSitesOnly-MaskedBadVariantsIndels/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz
outdir=$wd/genome-stats/pi/1MbWindowSize
mkdir -p $outdir
fai=~/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa.fai 
script=/u/home/a/ab08028/klohmueldata/annabel_data/genome-stats-scripts/pi/slidingWindowHeterozygosity.jar.ab.py
# run it:

python $script --vcf $vcf --window_size $window --step_size $step --fai $fai

# note that output is in: $wd/HQSitesOnly-MaskedBadVariantsIndels
# so move it to outdir
mv ${vcf}_het_${window}win_${step}step.txt $outdir

sleep 5m
