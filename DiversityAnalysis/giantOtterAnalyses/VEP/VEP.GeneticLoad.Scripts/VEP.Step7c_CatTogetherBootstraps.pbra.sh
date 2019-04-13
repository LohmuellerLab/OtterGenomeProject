# updated to work on new bootstraps (feb 2019) that use average that includes northern sea otter
###### 7c: cat together results 
date=20171206
spp=pbra

wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/bootstrapDerivedAlleles_avgInclNSO

# set up header:
head -n1 $wd/results/${spp}.group.1.10bootstraps.results.txt > $wd/${spp}.allBootstrapsConcatted.results.txt

for i in {1..100}
do
echo $i
rep=${spp}.group.${i}.10bootstraps.results.txt
grep -v "Missense" $wd/results/$rep >> $wd/${spp}.allBootstrapsConcatted.results.txt
# this grep -v "Missense" is just to avoid the header (all other lines only contain numbers)
# it has nothing to do with the missense category; that is just the first line of the header
done