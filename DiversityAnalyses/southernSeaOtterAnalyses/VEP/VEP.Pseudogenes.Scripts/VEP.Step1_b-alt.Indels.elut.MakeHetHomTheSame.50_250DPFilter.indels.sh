#!/bin/bash
#$ -l h_rt=100:00:00,h_data=2G,highp
#$ -N elutvepIndelsCDSOnly
#$ -pe shared 3
#$ -cwd
#$ -m bea
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -M ab08028
## Run VEP on homozygous and heterozygous snps:
##################### INDELs ########################
# to install:
# module load perl
# module load htslib/1.3.2
# run with --NO_HTSLIB
# say y to cache file; install 88 (mustela)
######### RUNNING ON CDS SITES ONLY ##########
# then load your modules:
source /u/local/Modules/default/init/modules.sh
module load perl/5.10.1
module load htslib
rundate=20171006 # date genotypes were called
PREFIX=01_Elut_CA_Gidget
vepdir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/ensembl-vep/
vcfdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter/cds_vcfs/indels/
# just running VEP on cds variants (hets/homs)
outdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter/vep-output/indels
mkdir -p $outdir
### heterozygotes:
echo "running VEP on heterozygotes"

hets=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.Heterozygous.vcf

### adding CANONICAL field so I can filter on that 
$vepdir/vep -v -i ${vcfdir}/$hets --fork 3 \
--cache --force_overwrite --species mustela_putorius_furo \
--numbers --domains --variant_class --canonical \
-o $outdir/${hets%.vcf}.VEP.output.minimalInfo.tbl
### homozygotes:
echo "running VEP on homozygotes"

homs=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.HomozygousAlt.vcf

$vepdir/vep -v -i ${vcfdir}/$homs --fork 3 \
--cache --force_overwrite --species mustela_putorius_furo \
--domains --numbers --variant_class --canonical \
-o ${outdir}/${homs%.vcf}.VEP.output.minimalInfo.tbl

sleep 5m
