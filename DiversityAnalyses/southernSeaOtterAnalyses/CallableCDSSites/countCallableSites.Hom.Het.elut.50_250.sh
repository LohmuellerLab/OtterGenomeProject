###### Count up callable sites in coding regions
date=20171006
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter
header=01_Elut_CA_Gidget
spp=elut
mkdir -p $wd/countCallableSites-CDSRegions
echo "chunk callableCDSSites" > $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.txt
echo "chunk callableCDSSites" > $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.VariantHets.txt
echo "chunk callableCDSSites" > $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.VariantHoms.txt
echo "chunk callableCDSSites" > $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.NonVariant.HomRefs.txt

for i in {1..323}
do
echo $i
vcf=$wd/cds_vcfs/${header}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf

tally=`grep -v "#" $vcf | sort | uniq | wc -l` # this will count up all non header lines in the cds vcf. Make sure that bedtools intersect
# does what you want it to. (it was actually adding duplicate lines, so now am adding sort | uniq to correct for that. I also changed the bedtools int script
# to use a merged bed file (0based) which should eliminate all duplicates. 
hets=`grep "0/1:" $vcf | sort | uniq | wc -l`
homs=`grep "1/1:" $vcf | sort | uniq | wc -l`
refs=`grep "0/0:" $vcf | sort | uniq | wc -l`
echo $i $tally >> $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.txt
echo $i $hets >> $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.VariantHets.txt
echo $i $homs >> $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.VariantHoms.txt
echo $i $refs >> $wd/countCallableSites-CDSRegions/${spp}.callableSites.ferretCDSRegions.NonVariant.HomRefs.txt

# unset variables so that they don't carry over
tally=
hets=
homs=
refs=
done

