## Final step get nucleotide genome coordinates of final functional ORs (or non functional ORs as well)
#
# functional list:
spp=$1

wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}
#cd $wd
#mkdir -p $wd/step_10_getFinalGenomeCoordsAndSeqs
funcIDs=$wd/step_9_FinalFunctionalList/${spp}.classI.and.classII.final.ORs.passedInspection.txt
finalNTSeqs=${spp}.step_10_results.finalORntSeqs.all.fna
# first step: pull out those protein seqs from step 5 fasta:
fasta_tool --select $funcIDs $finalNTSeqs > $wd/step_10_getFinalGenomeCoordsAndSeqs/${spp}.classI.and.classII.final.ORs.passedInspection.ntSeqs.fna
