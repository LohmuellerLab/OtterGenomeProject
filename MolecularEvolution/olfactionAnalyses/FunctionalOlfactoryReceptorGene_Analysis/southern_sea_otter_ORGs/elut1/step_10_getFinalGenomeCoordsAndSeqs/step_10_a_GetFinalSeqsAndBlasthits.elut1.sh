## Final step get nucleotide genome coordinates of final functional ORs (or non functional ORs as well)
#
# functional list:
spp=elut1

wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}
dupwd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}.dup
mkdir -p $wd/step_10_getFinalGenomeCoordsAndSeqs
db=$wd/step_1_ab_blastDb_blastResults/sea_otter_23May2016_bS9RH.fasta
 # this will be different for each spp
funcIDs=$wd/step_9_FinalFunctionalList/${spp}.classI.and.classII.final.ORs.passedInspection.txt
step5fasta=$wd/step_5_chooseM/step_5_result_mafftAlignment.pickedM.fasta
dupstep5fasta=${dupwd}/step_5_chooseM/step_5_result.pickedM.fasta

# need to also pull from dup ones.

# first step: pull out those protein seqs from step 5 fasta:
fasta_tool --select $funcIDs $step5fasta > $wd/step_10_getFinalGenomeCoordsAndSeqs/${spp}.classI.and.classII.final.ORs.passedInspection.AASeqs.faa
fasta_tool --select $funcIDs $dupstep5fasta >> $wd/step_10_getFinalGenomeCoordsAndSeqs/${spp}.classI.and.classII.final.ORs.passedInspection.AASeqs.faa
# doing it in two steps 
# use blast:
query=$wd/step_10_getFinalGenomeCoordsAndSeqs/${spp}.classI.and.classII.final.ORs.passedInspection.AASeqs.faa # all final protein sequences of functional ORs (after choose M in step 5)

# db has to be set for each script!

tblastn -query $query -db $db -num_threads 20 -outfmt 6 -max_target_seqs 1 > $wd/step_10_getFinalGenomeCoordsAndSeqs/step_10_finalFunctionalPGenes_blastOutput.txt
# Put Results in R

