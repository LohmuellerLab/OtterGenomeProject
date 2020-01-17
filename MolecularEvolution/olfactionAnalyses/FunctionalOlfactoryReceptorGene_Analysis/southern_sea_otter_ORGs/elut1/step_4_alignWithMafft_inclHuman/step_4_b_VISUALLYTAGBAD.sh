#use jalview to look for areas with greater than 5 gaps in each of the TM domains rlative to human sequence


#some will be so bad that they cause gaps in the ref sequence, get rid of them and realign 

# save names of sequences in a file :
badseqs=Step_4_badSeqsVisual_1.txt
rundate=20170925 # date you did alignment 
input=/work2/abeichman/gangLi_OlfactionPipeline/step_3_findORFs/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta
outputdir=/work2/abeichman/gangLi_OlfactionPipeline/step_4_alignWithMafft_inclHuman

# use maker's fasta_tool to remove them:
fasta_tool --remove $badseqs $input > $outputdir/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.REMOVEDbadSeqs1.fasta

