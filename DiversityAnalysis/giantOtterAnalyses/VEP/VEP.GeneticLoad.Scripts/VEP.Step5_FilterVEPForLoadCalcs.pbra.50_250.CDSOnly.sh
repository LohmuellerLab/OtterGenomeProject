#!/bin/bash
#$ -N filterPbraVEP
#$ -l h_data=5G,h_rt=10:00:00,highp
#$ -m bea
#$ -M ab08028
#$ -cwd
########## Want to filter Elut VEP Results to get the following:
### For pbra,  het/hom already separated 

# synonymous_variant 
# missense_variant (levels)
# homozygous/heterozygous

date=20171206 # date of vcfs
indir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/vep-output/
PREFIX=PteBra_1
# can pull out all sorts of things 
outdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/vep-output/filteredVepOutput2-ForLoadCalcs
mkdir -p $outdir
vepdir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/ensembl-vep/

#### Note: this only works for a tbl with a single individual! 
############## a. Separate hets/homs (already done for Pbra): ##########################

hetInput=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.Heterozygous.VEP.output.minimalInfo.tbl
homInput=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.VEP.output.minimalInfo.tbl

### Another addition: want to make sure it's in canonical protein so different sites dont' count twice.
# 20180206 : added CANONICAL is YES
############## b. Filter Synonymous/NS ##########################
###### NS (missense):
$vepdir/filter_vep --filter "Consequence is missense_variant and CANONICAL is YES" \
--input_file $indir/$homInput \
--output_file $outdir/${homInput%.tbl}.missense.tbl \
--force_overwrite
# note: only-matched only pulls the part of the line that matches so that if something 
# is annotated as >1 thing, only the matching part is included 

$vepdir/filter_vep --filter "Consequence is missense_variant and CANONICAL is YES" \
--input_file $indir/$hetInput \
--output_file $outdir/${hetInput%.tbl}.missense.tbl \
--force_overwrite

###### LOF (stop_gained):
$vepdir/filter_vep --filter "Consequence is stop_gained and CANONICAL is YES" \
--input_file $indir/$homInput \
--output_file $outdir/${homInput%.tbl}.stopgained.tbl \
--force_overwrite


$vepdir/filter_vep --filter "Consequence is stop_gained and CANONICAL is YES" \
--input_file $indir/$hetInput \
--output_file $outdir/${hetInput%.tbl}.stopgained.tbl \
--force_overwrite

###### LOF 2: stopgained in a DOMAIN
$vepdir/filter_vep --filter "Consequence is stop_gained and CANONICAL is YES and DOMAIN" \
--input_file $indir/$homInput \
--output_file $outdir/${homInput%.tbl}.stopgained.domain.tbl \
--force_overwrite


$vepdir/filter_vep --filter "Consequence is stop_gained and CANONICAL is YES and DOMAIN" \
--input_file $indir/$hetInput \
--output_file $outdir/${hetInput%.tbl}.stopgained.domain.tbl \
--force_overwrite

###### S (synonymous):

$vepdir/filter_vep --filter "Consequence is synonymous_variant and CANONICAL is YES" \
--input_file $indir/$homInput \
--output_file $outdir/${homInput%.tbl}.synonymous.tbl \
--force_overwrite

$vepdir/filter_vep --filter "Consequence is synonymous_variant and CANONICAL is YES" \
--input_file $indir/$hetInput \
--output_file $outdir/${hetInput%.tbl}.synonymous.tbl \
--force_overwrite
