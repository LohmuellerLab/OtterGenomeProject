# step a:
input1=step_2_results.olfacUniqueBlastHits.containNNN.stranded.fasta
input2=step_2_results.olfacUniqueBlastHits.within30ofEnds.stranded.fasta
input3=step_2_results.olfacUniqueBlastHits.notCloseToEnds.noNNN.stranded.fasta
perl step_3_a_Find_ORF-ABmod.generic.len120.pl $input1

perl step_3_a_Find_ORF-ABmod.generic.len120.pl $input2

perl step_3_a_Find_ORF-ABmod.generic.len120.pl $input3

# step b: 
perl step_3_b_Find_ORF-ABmod.WithPrintStatements.len120.generic.pl $input1
perl step_3_b_Find_ORF-ABmod.WithPrintStatements.len120.generic.pl $input2
perl step_3_b_Find_ORF-ABmod.WithPrintStatements.len120.generic.pl $input3



