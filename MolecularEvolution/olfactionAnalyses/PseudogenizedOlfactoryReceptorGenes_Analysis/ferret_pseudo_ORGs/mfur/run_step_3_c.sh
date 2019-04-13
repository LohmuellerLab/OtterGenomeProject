
input1=ORF_longThan_120_bp_step_2_results.olfacUniqueBlastHits.containNNN.stranded.fasta
input2=ORF_longThan_120_bp_step_2_results.olfacUniqueBlastHits.within30ofEnds.stranded.fasta
input3=ORF_longThan_120_bp_step_2_results.olfacUniqueBlastHits.notCloseToEnds.noNNN.stranded.fasta
sh step_3_c_GetAllORFIDs.Lengths.sh mfur $input1
sh step_3_c_GetAllORFIDs.Lengths.sh mfur $input2
sh step_3_c_GetAllORFIDs.Lengths.sh mfur $input3


# then in R get the longest of each site
