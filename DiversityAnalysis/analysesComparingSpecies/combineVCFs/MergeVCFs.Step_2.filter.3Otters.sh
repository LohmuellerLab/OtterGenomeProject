#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=20G
#$ -N mergeVCFs2
#$ -o /u/flashscratch/a/ab08028/otters/reports
#$ -e /u/flashscratch/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
#$ -t 1-323
# Combine variants 
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
# vcftools need to be able to find tabix; put it in path as well and then use -V to send that path to the job
rundate1=20171006 
rundate2=20171206
rundate3=20190215_nso_lib283only
PREFIX1=01_Elut_CA_Gidget
sample1=GIDGET
sample2=PteBra_1
sample3=ELFIN
PREFIX2=PteBra_1
PREFIX3=02_Elut_SEAK_Elfin
spp1="elut1_SSO"
spp2="pbra"
spp3="elut2_NSO"
# same ref for both (mfur)
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar

wd=/u/home/a/ab08028/kirk-bigdata/annabel/otters/
outdir=$wd/mergedVCFs_all3Spp/AllSites
snpdir=$wd/mergedVCFs_all3Spp/VariantsOnly

mkdir -p $outdir
mkdir -p $snpdir

# elut
dir1=/u/home/a/ab08028/kirk-bigdata/annabel/otters/vcfs/southern_sea_otter/HQSitesOnly-MaskedBadVariantsIndels/
vcf1=${PREFIX1}.raw_variants.${SGE_TASK_ID}.${rundate1}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz

# pbra:
dir2=/u/home/a/ab08028/kirk-bigdata/annabel/otters/vcfs/giant_otter/HQSitesOnly-MaskedBadVariantsIndels/
vcf2=${PREFIX2}.raw_variants.${SGE_TASK_ID}.${rundate2}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz

# elut2 (northern sea otter)
dir3=/u/flashscratch/a/ab08028/otters/vcfs/vcf_20190215_nso_lib283only_50_250_Filter/HQSitesOnly-MaskedBadVariantsIndels
vcf3=02_Elut_SEAK_Elfin.raw_variants.${SGE_TASK_ID}.20190215_nso_lib283only.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz

java -Xmx8g -jar $GATK \
-T SelectVariants \
-R $REFERENCE \
-V $outdir/combined_3Otters_${SGE_TASK_ID}.HQsites.Only.rmDotGenotypes.rmBadVars.CalledInAllOnly.bcftoolsMerged.vcf.gz \
--maxNOCALLfraction 0 \
--excludeNonVariants \
--restrictAllelesTo BIALLELIC \
-o $snpdir/combined_3Otters_${SGE_TASK_ID}.HQsites.Only.rmDotGenotypes.rmBadVars.CalledInAllOnly.bcftoolsMerged.VarsOnly.vcf.gz

sleep 5m
