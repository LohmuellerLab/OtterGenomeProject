#### run on SIRIUS:
spp=elut1.dup
# full genome including duplicate scaffs:
REFERENCE=/data3/abeichman/gidgetDeNovoGenome/DovetailAssembly/sea_otter_23May2016_bS9RH.fasta

# this doesn't have any length filter on hits, has +- 300 bp
bedPlus300=/work2/abeichman/gangLi_OlfactionPipeline_pgene_${spp}/step_1_c_CleanedByEvalue_ABScript_noLengthFilter/nonFunctionalORs.bed




# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 300 bp on either side (skipping all the stuff Gang did in Step 2 that made sure
# the query start and end matched the hit, instead just giving a big 200 bp buffer. Hopefully works
wd=/work2/abeichman/gangLi_OlfactionPipeline_pgene_${spp}/step_2_getSequencesFromGenome

cd $wd
bedtools getfasta -s -name -bed $bedPlus300 -fi $REFERENCE -fo step_2_results.olfacUniqueBlastHits.stranded.fasta
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 200bp to each side of your range
