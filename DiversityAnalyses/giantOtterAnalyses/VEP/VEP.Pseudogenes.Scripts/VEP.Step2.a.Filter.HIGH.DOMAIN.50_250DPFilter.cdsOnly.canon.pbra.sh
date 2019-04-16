#! /bin/bash
#$ -cwd
#$ -l h_rt=25:00:00,h_data=10G,highp
#$ -N filtervep
#$ -m abe
#$ -M ab08028

# Want to filter all results for HIGH impact (any kind) in a DOMAIN
rundate=20171206 # date of vcfs
vepdir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/ensembl-vep/
indir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter/vep-output/
outdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${rundate}_50_250DPFilter/vep-output/filteredVepOutput1_HighImpact
PREFIX=PteBra_1
mkdir -p $outdir
# want to filter SNPs, Indels, Homozygotes, Hets for HIGH and DOMAINS (skip NOT domain for now)
indelHet=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.Heterozygous.VEP.output.minimalInfo.tbl
indelHom=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.HomozygousAlt.VEP.output.minimalInfo.tbl

snpHet=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.Heterozygous.VEP.output.minimalInfo.tbl
snpHom=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.VEP.output.minimalInfo.tbl

# also adding: CANONICAL is YES
############## INDELS ####################
# filter out HIGH impact and DOMAINS again from final tables:
$vepdir/filter_vep --filter "IMPACT is HIGH and DOMAINS and CANONICAL is YES" \
--input_file $indir/indels/$indelHet \
--output_file $outdir/${indelHet%.tbl}.HIGH.DOMAIN.CANON.tbl \
--force_overwrite

$vepdir/filter_vep --filter "IMPACT is HIGH and DOMAINS and CANONICAL is YES" \
--input_file $indir/indels/$indelHom \
--output_file $outdir/${indelHom%.tbl}.HIGH.DOMAIN.CANON.tbl \
--force_overwrite

############## SNPs ####################
# filter out HIGH impact and DOMAINS again from final tables:
$vepdir/filter_vep --filter "IMPACT is HIGH and DOMAINS and CANONICAL is YES" \
--input_file $indir/$snpHet \
--output_file $outdir/${snpHet%.tbl}.HIGH.DOMAIN.CANON.tbl \
--force_overwrite

$vepdir/filter_vep --filter "IMPACT is HIGH and DOMAINS and CANONICAL is YES" \
--input_file $indir/$snpHom \
--output_file $outdir/${snpHom%.tbl}.HIGH.DOMAIN.CANON.tbl \
--force_overwrite

sleep 5m