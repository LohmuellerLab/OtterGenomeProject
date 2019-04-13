#!/bin/bash
#$ -N elut2Bootstraps
#$ -l h_data=20G,h_rt=10:00:00
#$ -m bea
#$ -M ab08028
#$ -cwd
#$ -t 1-100
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports

# doing 100 groups of 10
source /u/local/Modules/default/init/modules.sh
module load R
date=20190215_nso_lib283only

spp=elut2

wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/bootstrapDerivedAlleles

script=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/northern_sea_otter_mapping/TryWithJustOneLibrary_fromStep5On/IndelRealignment-BranchOfScripts/VEP_50_250DPFilter/includesSharedHomAltVars_UseForMainVEPResults/Steps5_7_geneticLoad_includesShared/VEP.Step7b_RunBootstrapsDerivedAlleles.${spp}_50_250DPFilter.R

Rscript $script ${SGE_TASK_ID}

sleep 5m