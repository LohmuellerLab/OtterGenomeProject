#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N abbaBaba
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -V
#$ -t 2-2
## run abba-baba
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts/test_abba_baba_3Otters

vcfdir=/u/home/a/ab08028/kirk-bigdata/annabel/otters/mergedVCFs_all3Spp/VariantsOnly
vcf=combined_3Otters_${SGE_TASK_ID}.HQsites.Only.rmDotGenotypes.rmBadVars.CalledInAllOnly.bcftoolsMerged.VarsOnly.vcf.gz
 # merged across all 3 otters (maybe per scaffold? or merged across all scaffolds; variant sites only)
outdir=/u/home/a/ab08028/kirk-bigdata/annabel/otters/abba-baba

python $scriptdir/abba.baba.AB.py --vcf $vcfdir/$vcf --treeOrder PteBra_1,GIDGET,Elfin --outPREFIX $SGE_TASK_ID --outdir $outdir
# can use --outPREFIX $SGE_TASK_ID if I do this split across scaffs. 
sleep 5m
