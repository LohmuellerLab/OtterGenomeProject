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
PREFIX=02_Elut_SEAK_Elfin
rundate=20190215_nso_lib283only # date of vcf folder
vcfdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}/
mkdir -p $vcfdir/stats
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/northern_sea_otter_mapping
# script dir okay
python $scriptdir/compute_DP.py --VCF $vcfdir/raw_variants/${PREFIX}.raw_variants.${SGE_TASK_ID}.${rundate}.vcf.gz --outfile $vcfdir/stats/DP.dist.scaff.${SGE_TASK_ID}.out

sleep 5m
