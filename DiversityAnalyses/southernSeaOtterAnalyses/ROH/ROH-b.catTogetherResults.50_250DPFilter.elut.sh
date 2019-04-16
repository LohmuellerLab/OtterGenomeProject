## step b:
# cat together all the ROHs:
date=20171006 # date of vcf calling
PREFIX=01_Elut_CA_Gidget
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/genome-stats/ROH
plinkoutdir=$wd/plinkOutputFiles/

mkdir -p $plinkoutdir/cattedResults
# set up header
head -n1 $plinkoutdir/${PREFIX}.raw_variants.1.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.Plink.out.hom > $plinkoutdir/cattedResults/${PREFIX}.scaffs1_220.catted.hom
for i in {1..220}
do
plinkHom=$plinkoutdir/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.Plink.out.hom
# grab chr name
grep -v "FID" $plinkHom >> $plinkoutdir/cattedResults/${PREFIX}.scaffs1_220.catted.hom
done 

# may also want to delete .summary files as they are huge! 
#rm $plinkoutdir/*summary 
