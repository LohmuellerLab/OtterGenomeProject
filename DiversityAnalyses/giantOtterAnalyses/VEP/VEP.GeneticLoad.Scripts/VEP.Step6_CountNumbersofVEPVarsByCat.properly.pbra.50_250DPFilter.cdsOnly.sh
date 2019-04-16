#### Get VEP counts 


date=20171206
PREFIX=PteBra_1
spp=pbra

indir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/vep-output/filteredVepOutput2-ForLoadCalcs
outputFile=${indir}/${PREFIX}.${date}.VepResultsCountSummary.txt
echo "species category genotype count" > $outputFile # this resets the output file every time 
for category in synonymous missense stopgained
do
for GT in HomozygousAlt Heterozygous
do
count="" # reset variable just in case 
# note addition of ".cdsSequence." below
inputFile=${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.${GT}.VEP.output.minimalInfo.${category}.tbl
count=`grep -v "#" $indir/$inputFile | awk '{print $2}' | sort | uniq | wc -l` # this excludes comment lines, gets the positions in format Chrom:Pos and makes sure they are a unique count
# you have do to the awk part because some variants affect mult. genes so the would be counted 2x. 
echo $spp $category $GT $count >> $outputFile
done
done
