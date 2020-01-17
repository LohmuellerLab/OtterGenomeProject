
# deduped ref genome:
REFERENCE=/data3/abeichman/gidgetDeNovoGenome/deduplicated99Genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta

bedPlus300=/work2/abeichman/gangLi_OlfactionPipeline/step_1_c_CleanedByEvalue_ABScript/cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed

# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)

wd=/work2/abeichman/gangLi_OlfactionPipeline/step_2_getSequencesFromGenome
cd $wd
bedtools getfasta -s -name -bed $bedPlus300 -fi $REFERENCE -fo step_2_results.olfacUniqueBlastHits.stranded.fasta
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 200bp to each side of your range 
