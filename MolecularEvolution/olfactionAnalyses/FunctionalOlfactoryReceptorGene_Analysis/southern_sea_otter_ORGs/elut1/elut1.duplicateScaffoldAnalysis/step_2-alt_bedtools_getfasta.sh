#### run on SIRIUS:

#don't use deduped ref genome here
# instead use full genome including duplicated scaffs (because getting those that were missed when genome was deduplicated)
REFERENCE=/data3/abeichman/gidgetDeNovoGenome/DovetailAssembly/sea_otter_23May2016_bS9RH.fasta

bedPlus300=/work2/abeichman/gangLi_OlfactionPipeline_elut1.dup/step_1_c_CleanedByEvalue_ABScript/elut1.DuplicateRemainderSeqs.cleanedByEvalue.250LengthFilter.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed


# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 300 bp on either side
wd=/work2/abeichman/gangLi_OlfactionPipeline_elut1.dup/step_2_getSequencesFromGenome
cd $wd
bedtools getfasta -s -name -bed $bedPlus300 -fi $REFERENCE -fo step_2_results.olfacUniqueBlastHits.stranded.fasta
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 300bp to each side of your range 
