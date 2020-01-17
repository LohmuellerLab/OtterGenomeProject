# look at alignment in jalview again and flag baddies
# put in : Step_4_badSeqsVisual_2.txt
badseqs=Step_4_badSeqsVisual_2.txt

# Realign again
input=/work2/abeichman/gangLi_OlfactionPipeline/step_4_alignWithMafft_inclHuman/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.REMOVEDbadSeqs1.fasta
outputdir=/work2/abeichman/gangLi_OlfactionPipeline/step_4_alignWithMafft_inclHuman

# use maker's fasta_tool to remove them:
fasta_tool --remove $badseqs $input > $outputdir/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.REMOVEDbadSeqs2.fasta

