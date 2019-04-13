#source /u/local/Modules/default/init/modules.sh
date=20171006
vcfdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/cds_vcfs/
PREFIX=01_Elut_CA_Gidget

### Only catting together variant files:
############## SNPs #######################
# get vcf header -- set up two files (homozygotes and heterozygotes)
# This is just grepping the "#" comment lines

grep "#" ${vcfdir}/${PREFIX}.raw_variants.1.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.vcf
grep "#" ${vcfdir}/${PREFIX}.raw_variants.1.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.Heterozygous.vcf
for i in {1..323}
do
echo $i
vcf=${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf

# get homozygotes:
grep "1/1" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.vcf
# get heterozygotes:
grep "0/1" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.Heterozygous.vcf
# get hom Ref: new addition 3/23/2018
grep "0/0" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.HomozygousReference.vcf

done

# check that it adds up

