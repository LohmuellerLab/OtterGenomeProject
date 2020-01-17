# after you've visually gotten rid of bad sequences, realign with mafft to see if you need to remove any more:
spp=$1

input=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_4_alignWithMafft_inclHuman/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.REMOVEDbadSeqs1.fasta
outputdir=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_4_alignWithMafft_inclHuman
rundate=`date +%Y%m%d`
# installed mafft v. 7.310
echo "mafft parameters: mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder" > $outputdir/mafftParams.${rundate}.txt
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $input > $outputdir/step_4_result_mafftAligment.${rundate}.removedBadSeqs1.fasta
