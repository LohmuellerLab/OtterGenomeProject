#!/bin/bash
#$ -N elutBootstraps
#$ -l h_data=20G,h_rt=10:00:00,highp
#$ -m bea
#$ -M ab08028
#$ -cwd
#$ -t 1-100
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
# doing 100 groups of 10
source /u/local/Modules/default/init/modules.sh
module load R
date=20171006
spp=elut

wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/bootstrapDerivedAlleles

script=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts/redo_elut1pbra_GeneticLoadBootstraps_avginclNSO/VEP.Step7b_RunBootstrapsDerivedAlleles.${spp}_50_250DPFilter.avgIncl.NSO.R

Rscript $script ${SGE_TASK_ID}

sleep 5m