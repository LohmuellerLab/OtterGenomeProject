#! /bin/bash
#$ -cwd
#$ -l h_rt=3:00:00,h_data=2G
#$ -N msmc
#$ -pe shared 16
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028


# run MSMC with default settings (for now)

source /u/local/Modules/default/init/modules.sh
module load python/3.4
msmc_tools=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc-tools
msmc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc
date=20190215_nso_lib283only # vcf date
rundate=`date +%Y%m%d` # msmc rundate 
spp=elut2
OUTDIR=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/msmcAnalysis/output_${rundate}
mkdir -p $OUTDIR

INPUTDIR=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/msmcAnalysis/inputFiles

$msmc -t 16 -o $OUTDIR/${spp}.msmc.out $INPUTDIR/chunk_*_postMultiHetSep.txt

