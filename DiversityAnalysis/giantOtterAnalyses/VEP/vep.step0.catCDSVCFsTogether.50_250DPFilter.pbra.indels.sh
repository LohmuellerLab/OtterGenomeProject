#source /u/local/Modules/default/init/modules.sh
date=20171206
vcfdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/cds_vcfs/indels
PREFIX=PteBra_1

### Only catting together variant files:
############## INDELS #######################
# get vcf header -- set up two files (homozygotes and heterozygotes)
# This is just grepping the "#" comment lines

grep "#" ${vcfdir}/${PREFIX}.raw_variants.1.${date}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.CDSOnly.fromMergedBed.vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.vcf

grep "#" ${vcfdir}/${PREFIX}.raw_variants.1.${date}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.CDSOnly.fromMergedBed.vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.Heterozygous.vcf

for i in {1..323}
do
echo $i
vcf=${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.CDSOnly.fromMergedBed.vcf
# get homozygotes:
grep "1/1" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.HomozygousAlt.vcf
# get heterozygotes:
grep "0/1" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.Heterozygous.vcf
done
# check that it adds up

