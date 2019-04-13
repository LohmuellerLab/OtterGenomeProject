#Prepare files for Pseudopipe on SIRIUS 
#What you need:
#repeat masked genome (done): giantOtterAnnotation_Feb2017/repeatMaskGenomeFasta_20170307/PteBra.a.lines.15kbOnly.REPEATMASKED.NNs.fasta
#unmasked genome with 1 file per scaffold: /work2/abeichman/giantOtterAnnotation_Feb2017/prepForPseudoPipe_20170307/chunkedGiantOtter15kb
#1 file per scaffold giving exon st/stop in 3rd/4th columns
#peptides : 
#giantOtterAnnotation_Feb2017/PteBra.a.lines.maker.output/makerRound1Output/finalRound1Files/final_PteBra.a.lines.all.maker.proteins.round1.renamed.fasta
#*** using 15kb fasta for elut and pbra because that was what was annotated! *** 
# 1. take genome with >15kb (what was annotated) and split using maker's fasta_tool into 1267 files (elut); 100 groups of 260 files (pbra)
#pbra: chunk first, move to sub dirs, then split (this is better than the manual way I did elut)
mkdir dna
cd dna
fasta_tool --chunks 100 PteBra.a.lines.15kbOnly.fasta # this splits into 100 chunks (000 ... 099)
for i in {000..099}
do
echo $i
mkdir group_${i}/
mv PteBra.a.lines.15kbOnly_${i}.fasta group_${i}/
cd group_${i}/
fasta_tool --split PteBra.a.lines.15kbOnly_${i}.fasta
cd ../
done
# then got rid of chunks, just keeping per scaffold splits: rm group_0*/PteBra.a.*
# now each group just has "flattened_line whatever" files in it, one per scaffold.
# 2. Take gff and extract exons; make start and stop pos 3rd and 4th columns
### pbra  ### use this in downstream steps; use pbra in same downstream steps (did these in their own directories)
finalgff=/work2/abeichman/giantOtterAnnotation_Feb2017/PteBra.a.lines.maker.output/makerRound1Output/finalRound1Files/final_PteBra.a.lines.all.round1.renamed.gff
## pull out exons
awk '{OFS="\t";if($3=="exon") print $1,$3,$4,$5,$6,$7,$8,$9}' $finalgff > allScaffolds.exons.CorrectColumns.pbra.txt


# 3. split up exons by scaffold. 

mkdir mysql 
for i in {000..099}
do
echo $i
scaffs=`ls dna/group_${i} | sed 's/.fasta/ /g' | tr '\n' ' '` # fixed this 20170320
mkdir mysql/group_${i}
for j in $scaffs
do
echo $j
awk -v scaff=$j '{if($1==scaff) print}' allScaffolds.exons.CorrectColumns.pbra.txt > mysql/group_${i}/chr${j}_exLocs
done
done
#  transfer it all to hoffman:
