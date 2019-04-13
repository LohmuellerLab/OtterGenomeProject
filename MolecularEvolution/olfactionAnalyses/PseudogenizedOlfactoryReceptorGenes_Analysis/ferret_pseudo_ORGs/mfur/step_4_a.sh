# cat together human ORJ sequence with the putatively truncated genes 
human=Human_OR2J3.fasta

cat $human ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.containNNN.stranded.LONGEST.fasta \
ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.within30ofEnds.stranded.LONGEST.fasta \
> putativelyTruncated.addedHuman.fasta


