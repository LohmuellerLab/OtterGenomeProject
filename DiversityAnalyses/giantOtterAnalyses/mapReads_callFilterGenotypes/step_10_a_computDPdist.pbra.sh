#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=10G,highp
#$ -N computeDP
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-50
# get DP dist for a sampling of top 50 scaffolds
source /u/local/Modules/default/init/modules.sh
module load python/2.7
PREFIX=PteBra_1
date=20171206 # date of vcf folder
vcfdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}
mkdir -p $vcfdir/stats
python compute_DP.py --VCF $vcfdir/${PREFIX}.raw_variants.${SGE_TASK_ID}.${date}.vcf.gz --outfile $vcfdir/stats/DP.dist.scaff.${SGE_TASK_ID}.out
sleep 5m
