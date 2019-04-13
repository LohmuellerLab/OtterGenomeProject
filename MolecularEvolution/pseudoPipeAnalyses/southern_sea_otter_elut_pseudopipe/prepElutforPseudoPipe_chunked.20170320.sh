# 1. take genome with >15kb (what was annotated) and chunk using maker's fasta_tool into 10 chunks:
# CHUNK BOTH REPEAT MASKED and NOT repeatmasked versions
# not masked:
fasta_tool --chunks 10 sea_otter_23May2016_bS9RH.15kbONLY.fasta 
mkdir unmaskedChunks 
mv sea_otter_23May2016_bS9RH.15kbONLY_*.fasta unmaskedChunks/
## NOTE!!! these chunks are just for splitting. You don't want to use the chunks for the reference in ppipe
# instead you need equivalent chunks from the *repeat masked* genome! Slightly confusing.
# Individual scaffold fastas should be not repeat masked; each chunk should be RM masked. So have to chunk twice.
# masked:
fasta_tool --chunks 10 sea_otter_23May2016_bS9RH.15kbONLY.RepeatRunner.RepeatMasker.MASKED.fasta 
mkdir maskedChunks 
mv sea_otter_23May2016_bS9RH.15kbONLY.RepeatRunner.RepeatMasker.MASKED*.fasta maskedChunks/

mkdir dna

for i in {00..09}
do
echo $i
mkdir dna/group_${i}/
# get masked chunk in there:
cp maskedChunks/sea_otter_23May2016_bS9RH.15kbONLY.RepeatRunner.RepeatMasker.MASKED_${i}.fasta dna/group_${i}/
# get unmasked chunk in there: 
cp unmaskedChunks/sea_otter_23May2016_bS9RH.15kbONLY_${i}.fasta dna/group_${i}/
cd dna/group_${i}/
# split **UNMASKED CHUNK *** 
fasta_tool --split sea_otter_23May2016_bS9RH.15kbONLY_${i}.fasta

cd ../../
done
# but then don't want to keep sea_otter_23May2016_bS9RH.15kbONLY_${i}.fasta in there, so remove it:
ls dna/group_*/sea_otter_23May2016_bS9RH.15kbONLY_*.fasta
# REMOVE: rm dna/group_*/sea_otter_23May2016_bS9RH.15kbONLY_*.fasta

# 2. Take gff and extract exons; make start and stop pos 3rd and 4th columsn (already done)
finalgff=/data3/abeichman/AnnotationAugust2016/round5MakerOutput_AED1.0_useThis/finalMakerOutput/final_sea_otter_23May2016_bS9RH.round5.AED1.0.blast1e-06.renamed.gff
awk '{OFS="\t";if($3=="exon") print $1,$3,$4,$5,$6,$7,$8,$9}' $finalgff > allScaffolds.exons.CorrectColumns.elut.txt
# 3. split up exons by scaffold. 
mkdir mysql 
for i in {00..09}
do
echo $i
scaffs=`ls dna/group_${i}| sed 's/.fasta/ /g' | tr '\n' ' '` # fixed this 20170320
# this generates an extra empty exLocs. 
mkdir mysql/group_${i}
for j in $scaffs
do
echo $j
awk -v scaff=$j '{if($1==scaff) print}' allScaffolds.exon.CorrectColumns.elut.txt > mysql/group_${i}/chr${j}_exLocs
done
done
# get rid of extra stupid ex locs made from genome chunks. 
rm mysql/group_0*/chrsea_otter*
# 4. transfer it all to hoffman: 
# unmasked dna files (1 per scaffold)
rsync -rz -v --ignore-existing --partial /work2/abeichman/prepElut_Forpseudopipe_chunked/dna ab08028@hoffman2.idre.ucla.edu:/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes/ppipe_input/enhydra_lutris_chunked
# exon location files:
rsync -rz -v --ignore-existing --partial /work2/abeichman/prepElut_Forpseudopipe_chunked/mysql ab08028@hoffman2.idre.ucla.edu:/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes/ppipe_input/enhydra_lutris_chunked

