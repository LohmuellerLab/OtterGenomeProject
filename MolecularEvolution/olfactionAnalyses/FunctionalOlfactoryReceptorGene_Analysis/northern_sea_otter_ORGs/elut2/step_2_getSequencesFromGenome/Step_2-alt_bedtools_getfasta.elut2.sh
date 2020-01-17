#### run on SIRIUS:

# deduped ref genome:
spp="elut2"
REFERENCE=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_1_ab_blastDb_blastResults/GCA_002288905.2_ASM228890v2_genomic.fna



bedPlus300=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_1_c_CleanedByEvalue_ABScript/${spp}.cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed

# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 300 bp on either side
wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_2_getSequencesFromGenome
cd $wd
bedtools getfasta -s -name -bed $bedPlus300 -fi $REFERENCE -fo step_2_results.olfacUniqueBlastHits.stranded.fasta
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 300bp to each side of your range 
