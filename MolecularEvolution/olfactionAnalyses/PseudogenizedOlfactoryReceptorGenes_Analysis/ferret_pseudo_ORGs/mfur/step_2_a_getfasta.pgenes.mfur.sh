#### run on SIRIUS:
spp=mfur
# deduped ref genome:
REFERENCE=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_1_ab_blastDb_blastResults/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa
# this doesn't have any length filter on hits, has +- 300 bp
bedPlus300=/work2/abeichman/gangLi_OlfactionPipeline_pgene_${spp}/step_1_c_CleanedByEvalue_ABScript_noLengthFilter/nonFunctionalORs.${spp}.bed




# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 300 bp on either side
wd=/work2/abeichman/gangLi_OlfactionPipeline_pgene_${spp}/step_2_getSequencesFromGenome

cd $wd
bedtools getfasta -s -name -bed $bedPlus300 -fi $REFERENCE -fo step_2_results.olfacUniqueBlastHits.stranded.fasta
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 300bp to each side of your range
