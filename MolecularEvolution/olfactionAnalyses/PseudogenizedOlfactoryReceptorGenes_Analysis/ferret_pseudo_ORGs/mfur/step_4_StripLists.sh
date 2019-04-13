# strip script:
files=`ls`


for file in $files
do
grep -v "#" $file | sed 's/Annotation for //g' | sed -E 's/-[0-9]+_aa.*$//g' > ${file%.txt}.IDList.format.txt
done

