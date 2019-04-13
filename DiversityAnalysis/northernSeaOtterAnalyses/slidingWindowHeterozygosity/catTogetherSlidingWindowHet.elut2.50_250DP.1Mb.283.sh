window=1000000
dirlabel=1MbWindowSize 
step=20000
date=20190215_nso_lib283only
 # date of vcf calling
PREFIX=02_Elut_SEAK_Elfin
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/genome-stats/pi

mkdir -p $wd/scaff1-220_catted

cat $wd/$dirlabel/${PREFIX}.raw_variants.1.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz_het_${window}win_${step}step.txt \
> $wd/scaff1-220_catted/${PREFIX}.raw_variants.1-220Scaffs.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz_het_${window}win_${step}step.txt
 
 
for i in {2..220}
do
echo $i
grep -v "chromo" $wd/$dirlabel/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz_het_${window}win_${step}step.txt \
>> $wd/scaff1-220_catted/${PREFIX}.raw_variants.1-220Scaffs.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz_het_${window}win_${step}step.txt
done
