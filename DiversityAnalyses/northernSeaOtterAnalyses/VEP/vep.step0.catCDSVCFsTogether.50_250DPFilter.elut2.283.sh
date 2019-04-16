#source /u/local/Modules/default/init/modules.sh
rundate=20190215_nso_lib283only
vcfdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/cds_vcfs/
PREFIX=02_Elut_SEAK_Elfin

### Only catting together variant files:
############## SNPs #######################
# get vcf header -- set up two files (homozygotes and heterozygotes)
# This is just grepping the "#" comment lines

vcf=${PREFIX}.raw_variants.1.${rundate}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf


# set up files with new headers 
grep "#" $vcfdir/$vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.vcf
grep "#" $vcfdir/$vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.Heterozygous.vcf
grep "#" $vcfdir/$vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.HomozygousReference.vcf

for i in {1..323}
do
echo $i
vcf=${PREFIX}.raw_variants.${i}.${rundate}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf

# get homozygotes:
grep "1/1" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.vcf
# get heterozygotes:
grep "0/1" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.Heterozygous.vcf
# get hom Ref: 
grep "0/0" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${rundate}.HQsites.Only.HomozygousReference.vcf
done

# check that it adds up
