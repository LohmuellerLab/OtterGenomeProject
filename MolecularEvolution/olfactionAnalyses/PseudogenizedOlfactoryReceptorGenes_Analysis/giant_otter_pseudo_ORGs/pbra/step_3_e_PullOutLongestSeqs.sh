# 3e

input1=ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.containNNN.stranded
input2=ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.within30ofEnds.stranded
input3=ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.notCloseToEnds.noNNN.stranded

fasta_tool --select ${input1}.IDlist.LONGEST ./${input1}.fasta > ${input1}.LONGEST.fasta
fasta_tool --select ${input2}.IDlist.LONGEST ./${input2}.fasta > ${input2}.LONGEST.fasta
fasta_tool --select ${input3}.IDlist.LONGEST ./${input3}.fasta > ${input3}.LONGEST.fasta
