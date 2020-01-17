# step 4 b: after visually tagging bad sequences, put names in a file (4a_Exclude.txt)
# use maker's fasta_tool to remove them:
# save names of sequences in a file :
# usage: script species dateOfAlignment
# e.g. script mfur 20171009
# e.g. script pbra 20171009
badseqs=4a_Exclude.txt
rundate=$2 # date you did alignment 
spp=$1
input=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_3_findORFs/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.fasta
outputdir=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_4_alignWithMafft_inclHuman

# use maker's fasta_tool to remove them:
fasta_tool --remove $badseqs $input > $outputdir/ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.addedHuman_0R2J3.REMOVEDbadSeqs1.fasta

