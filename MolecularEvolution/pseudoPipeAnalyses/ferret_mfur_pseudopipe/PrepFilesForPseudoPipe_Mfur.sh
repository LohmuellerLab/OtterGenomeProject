#Prepare files for Pseudopipe on SIRIUS 
#What you need:
#repeat masked genome (done): 
#unmasked genome with 1 file per scaffold: 
#1 file per scaffold giving exon st/stop in 3rd/4th columns
#peptides : Mustela_putorius_furo.MusPutFur1.0.pep.all.fa

#
# doing the same for mfur for consistency.
# 1. take genome with >15kb (what was annotated) and split using maker's fasta_tool

###  MFUR: 
# pwd: /work2/abeichman/prepMFUR_forPpipe/ with repeat masked and regular 15kb genomes in there
mkdir dna
cd dna
fasta_tool --chunks 100 ../Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.15kb.fa # this splits into 100 chunks (000 ... 099)
# ended up split into 83 groups 
for i in {000..083}
do
echo $i
mkdir group_${i}/
mv ../Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.15kb_${i}.fasta group_${i}/
cd group_${i}/
fasta_tool --split Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.15kb_${i}.fasta
cd ../
done
# then got rid of chunks, just keeping per scaffold splits: rm group_0*/PteBra.a.*
# now each group just has "flattened_line whatever" files in it, one per scaffold. Nice! 
# 2. Take gff and extract exons; make start and stop pos 3rd and 4th columns
### mfur  ### use this in downstream steps
finalgff=/work2/abeichman/prepMFUR_forPpipe/Mustela_putorius_furo.MusPutFur1.0.89.gff3

# 20170717: downloaded from ftp://ftp.ensembl.org/pub/release-89/gff3/mustela_putorius_furo/Mustela_putorius_furo.MusPutFur1.0.89.gff3.gz
## pull out exons
awk '{OFS="\t";if($3=="exon") print $1,$3,$4,$5,$6,$7,$8,$9}' $finalgff > allScaffolds.exons.CorrectColumns.mfur.txt


# 3. split up exons by scaffold. 
# so instead, get scaff names for each group.
# pwd: /work2/abeichman/prepMFUR_forPpipe/ with repeat masked and regular 15kb genomes in there
mkdir mysql 
for i in {000..083}
do
echo $i
scaffs=`ls dna/group_${i} | sed 's/.fasta/ /g' | tr '\n' ' '` # fixed this 20170320
mkdir mysql/group_${i}
for j in $scaffs
do
echo $j
awk -v scaff=$j '{if($1==scaff) print}' allScaffolds.exons.CorrectColumns.mfur.txt > mysql/group_${i}/chr${j}_exLocs
done
done
# 5. repeat masked genome -- need section of it in each group for blasting

#repeat masked genome: /work2/abeichman/prepMFUR_forPpipe/Mustela_putorius_furo.MusPutFur1.0.dna_rm.toplevel.RepeatMaskedByEnsembl.fa
#downloaded from: ftp://ftp.ensembl.org/pub/release-89/fasta/mustela_putorius_furo/dna/ on july 7, 2017
# select only 15kb scaffolds on RM genome:
./filter_fasta_by_seq_length.pl Mustela_putorius_furo.MusPutFur1.0.dna_rm.toplevel.RepeatMaskedByEnsembl.fa Mustela_putorius_furo.MusPutFur1.0.dna_rm.toplevel.RepeatMaskedByEnsembl.15kbOnly.fa  15000
# then do the same --chunks approach (100 ends up with 83)
fasta_tool --chunks 100 Mustela_putorius_furo.MusPutFur1.0.dna_rm.toplevel.RepeatMaskedByEnsembl.15kbOnly.fa
# then sort them into the dna groups:
for i in {000..083}
do
mv Mustela_putorius_furo.MusPutFur1.0.dna_rm.toplevel.RepeatMaskedByEnsembl.15kbOnly_${i}.fasta dna/group_${i}
done


### need to have Repeat Masked genome in each group -- don't do it until you're on Hoffman
# transfer it all to hoffman: 

# unmasked dna files (1 per scaffold)
rsync -rz -v --ignore-existing --partial /work2/abeichman/prepMFUR_forPpipe/dna ab08028@hoffman2.idre.ucla.edu:/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes/ppipe_input/mustela_putorius_furo
# exon location files:
rsync -rz -v --ignore-existing --partial /work2/abeichman/prepMFUR_forPpipe/mysql ab08028@hoffman2.idre.ucla.edu:/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes/ppipe_input/mustela_putorius_furo
#
