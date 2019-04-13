# want to find sequences that are within 300 bp of a gap (NNs)
# could set some limit of NNs, but I find that even two are indicative of the presence of a gap in
# the otter assembly. Say 3 NNNs consecutively?
spp=$1
wd=/work2/abeichman/gangLi_OlfactionPipeline_pgene_${spp}
cd $wd/step_2_getSequencesFromGenome
mkdir -p closeToGaps
# use fasta_tool
fasta_tool --grep_seq "NNN" step_2_results.olfacUniqueBlastHits.stranded.fasta > closeToGaps/step_2_results.olfacUniqueBlastHits.containNNN.stranded.fasta

# this pulls sequences that have an NNN somewhere in the hit or within 300bp of it

mkdir -p closeToEnds 
# this will pull out sequences close to end of scaffolds (Niimura wants end of contig though)
# in step 1c I flagged the seqs within 30bp of end of scaffold:
within30ofEnds=$wd/step_1_c_CleanedByEvalue_ABScript_noLengthFilter/regionsCloseToStartEndofScaffold.txt

fasta_tool --select $within30ofEnds  step_2_results.olfacUniqueBlastHits.stranded.fasta > closeToEnds/step_2_results.olfacUniqueBlastHits.within30ofEnds.stranded.fasta

# left overs: ("other")
mkdir -p other
fasta_tool --grepv_seq "NNN" step_2_results.olfacUniqueBlastHits.stranded.fasta | fasta_tool --remove $within30ofEnds /dev/stdin > other/step_2_results.olfacUniqueBlastHits.notCloseToEnds.noNNN.stranded.fasta
