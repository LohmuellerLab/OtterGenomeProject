spp=elut #not elut1 bc grepping in fasta file
input1=ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.containNNN.stranded.fasta
input2=ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.within30ofEnds.stranded.fasta
input3=ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.notCloseToEnds.noNNN.stranded.fasta
sh step_3_c_GetAllORFIDs.Lengths.sh $spp $input1
sh step_3_c_GetAllORFIDs.Lengths.sh $spp $input2
sh step_3_c_GetAllORFIDs.Lengths.sh $spp $input3


# then in R get the longest of each site
