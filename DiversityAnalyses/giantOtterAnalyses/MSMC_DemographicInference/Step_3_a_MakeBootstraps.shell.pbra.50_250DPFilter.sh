### Want to figure out bootstrapping MSMC
module load python/3.4
# this is really fast, can do it in the shell 
# do 20 bootstraps:

PREFIX=PteBra_1
rundate=20171206
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter
msmc_tools=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc-tools
msmc=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc

# chunk size: default 5000000
# chunks per chrom: how many chunks on a boot strap chrom (20)
# number of chromosomes (30)
# n: num of bootstraps (20)
# seed: initialize random seed  
cd $wd/msmcAnalysis
mkdir -p $wd/msmcAnalysis/bootstraps
# doing 20 chromosomes, (each of 100Mb) so it's 2Gb of sequence, same as the 220 scaffolds
$msmc_tools/multihetsep_bootstrap.py -n 20 --nr_chromosomes 20 bootstraps/bootstrap-input/bootstrap_number $wd/msmcAnalysis/inputFiles/*txt
