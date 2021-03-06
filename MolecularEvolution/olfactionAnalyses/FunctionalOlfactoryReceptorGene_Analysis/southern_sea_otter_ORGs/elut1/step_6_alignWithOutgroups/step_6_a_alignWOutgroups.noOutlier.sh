# you picked the best starting M for each sequence
# then use step 5b to get rid of an extra human sequence, and add back in original human sequence
# and the outgroups. 
# then align and make a tree 

input=/work2/abeichman/gangLi_OlfactionPipeline/step_5_chooseM/step_5_result_mafftAlignment.pickedM.wHuman.wOutgroups.wRepresentatives.noOutlier.fasta

outputdir=/work2/abeichman/gangLi_OlfactionPipeline/step_6_alignWithOutgroups

rundate=`date +%Y%m%d`
# installed mafft v. 7.310
echo "mafft parameters: mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder" > $outputdir/mafftParams.${rundate}.txt
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $input > $outputdir/step_6_result_mafftAligment.wOutgroups.wRepresentatives.noOutlier.${rundate}.fasta

