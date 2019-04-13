#### Get VEP counts 


date=20190215_nso_lib283only # date of vcf calling
PREFIX=02_Elut_SEAK_Elfin
spp=elut2

indir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250_Filter/vep-output/filteredVepOutput2-ForLoadCalcs
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
# you have do to the awk part because some variants affect mult. genes so would be counted 2x. 
echo $spp $category $GT $count >> $outputFile
done
done
