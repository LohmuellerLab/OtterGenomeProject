## step 4:

# going to align with mafft follow gang li's recs:

spp=$1
# added the human OR2J3.fas sequence to the results of step 3 manually
inputdir=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_3_findORFs

cat $inputdir/Human_OR2J3.fasta $inputdir/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta > $inputdir/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta

mafftInput=$inputdir/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta

outputdir=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_4_alignWithMafft_inclHuman
rundate=`date +%Y%m%d`
# installed mafft v. 7.310
echo "mafft parameters: mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder" > $outputdir/mafftParams.${rundate}.txt
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $mafftInput > $outputdir/step_4_result_mafftAligment.${rundate}.fasta