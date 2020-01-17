### after step 7, you plotted in R, pulled out any weird seqs and wrote out lists of the class I and classI ORs
# want to align them separately, with the representative sequences:
wd=/work2/abeichman/gangLi_OlfactionPipeline/step_8_alignClassesSeparately
classIlist=$wd/classI.tips.includingReps.txt
classIIlist=$wd/classII.tips.includingReps.txt
inputfasta=/work2/abeichman/gangLi_OlfactionPipeline/step_5_chooseM/step_5_result_mafftAlignment.pickedM.wHuman.wOutgroups.wRepresentatives.noOutlier.fasta

fasta_tool --select $classIlist $inputfasta > $wd/classI.wReps.fasta
fasta_tool --select $classIIlist $inputfasta > $wd/classII.wReps.fasta 

rundate=`date +%Y%m%d`
# installed mafft v. 7.310
echo "mafft parameters: mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder" > $wd/mafftParams.${rundate}.txt
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $wd/classI.wReps.fasta  > $wd/step_8_result_mafftAligment.classI.${rundate}.fasta
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $wd/classII.wReps.fasta  > $wd/step_8_result_mafftAligment.classII.${rundate}.fasta
