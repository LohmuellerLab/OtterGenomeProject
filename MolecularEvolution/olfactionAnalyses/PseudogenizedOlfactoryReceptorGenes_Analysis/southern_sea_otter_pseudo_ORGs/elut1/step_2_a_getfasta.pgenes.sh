#### run on SIRIUS:
spp=elut1
# deduped ref genome:
REFERENCE=/data3/abeichman/gidgetDeNovoGenome/deduplicated99Genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta

# this doesn't have any length filter on hits, has +- 300 bp
bedPlus300=/work2/abeichman/gangLi_OlfactionPipeline_pgene_elut1/step_1_c_CleanedByEvalue_ABScript_noLengthFilter/nonFunctionalORs.bed




# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 300 bp on either side
wd=/work2/abeichman/gangLi_OlfactionPipeline_pgene_${spp}/step_2_getSequencesFromGenome

cd $wd
bedtools getfasta -s -name -bed $bedPlus300 -fi $REFERENCE -fo step_2_results.olfacUniqueBlastHits.stranded.fasta
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 300bp to each side of your range
