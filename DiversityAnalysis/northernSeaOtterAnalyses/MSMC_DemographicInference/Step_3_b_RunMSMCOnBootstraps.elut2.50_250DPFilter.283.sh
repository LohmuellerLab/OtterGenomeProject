#! /bin/bash
#$ -cwd
#$ -l h_rt=3:00:00,h_data=2G,highp
#$ -N bootstrapMSMC
#$ -pe shared 10
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-20
source /u/local/Modules/default/init/modules.sh

module load python/3.4

PREFIX=02_Elut_SEAK_Elfin
rundate=20190215_nso_lib283only
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250_Filter
msmc_tools=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc-tools
msmc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc
mkdir -p $wd/msmcAnalysis/bootstraps/bootstrap-output
# array of the number of bootstraps
$msmc -t 10 -o $wd/msmcAnalysis/bootstraps/bootstrap-output/bootstrap_number_${SGE_TASK_ID}.out $wd/msmcAnalysis/bootstraps/bootstrap-input/bootstrap_number_${SGE_TASK_ID}/*txt

sleep 5m
