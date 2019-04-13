## step 4:
# added the human OR2J3.fas sequence to the results of step 3 manually
# going to align with mafft follow gang li's recs:

spp=elut1.dup
wd=/work2/abeichman/gangLi_OlfactionPipeline_pgene_${spp}/step_4_align

input=$wd/putativelyTruncated.addedHuman.fasta
outputdir=$wd

rundate=`date +%Y%m%d`
# installed mafft v. 7.310
echo "mafft parameters: mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder" > outputdir/mafftParams.${rundate}.txt
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $input > $outputdir/step_4_result_mafftAligment.${rundate}.putTruncated.fasta
