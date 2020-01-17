### after step 7, you plotted in R, pulled out any weird seqs and wrote out lists of the class I and classI ORs
# want to align them separately, with the representative sequences:
spp=$1
wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_8_alignClassesSeparately
classIlist=$wd/classI.tips.includingReps.noOutlier.txt
classIIlist=$wd/classII.tips.includingReps.noOutlier.txt
inputfasta=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_5_chooseM/step_5_result.pickedM.wHuman.wOutgroups.wRepresentatives.fasta # don't need to remove outlier from this file
# because it isn't included in the class lists 

fasta_tool --select $classIlist $inputfasta > $wd/classI.wReps.fasta
fasta_tool --select $classIIlist $inputfasta > $wd/classII.wReps.fasta 

rundate=`date +%Y%m%d`
# installed mafft v. 7.310
echo "mafft parameters: mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder" > $wd/mafftParams.${rundate}.txt
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $wd/classI.wReps.fasta  > $wd/step_8_result_mafftAligment.classI.${rundate}.fasta
mafft --preservecase --genafpair --maxiterate 16 --thread 16 --inputorder $wd/classII.wReps.fasta  > $wd/step_8_result_mafftAligment.classII.${rundate}.fasta
