# concat ABBA-BABA results 

wd=/u/home/a/ab08028/kirk-bigdata/annabel/otters/abba-baba
mkdir -p $wd/concattedResults
# get header
grep scaffold $wd/ABBA-BABA.1.txt > $wd/concattedResults/ABBA-BABA.allScaffs.txt
for i in {1..323}
do
echo $i
# want to skip the header on each one
grep -v scaffold $wd/ABBA-BABA.${i}.txt >> $wd/concattedResults/ABBA-BABA.allScaffs.txt
done

gzip -f $wd/concattedResults/ABBA-BABA.allScaffs.txt