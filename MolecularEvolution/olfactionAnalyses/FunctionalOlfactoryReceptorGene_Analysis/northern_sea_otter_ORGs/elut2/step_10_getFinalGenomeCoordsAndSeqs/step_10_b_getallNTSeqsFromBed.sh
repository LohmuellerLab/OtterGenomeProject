#### run on SIRIUS:

# deduped ref genome:
spp="elut2"
wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}

REFERENCE=$wd/step_1_ab_blastDb_blastResults/GCA_002288905.2_ASM228890v2_genomic.fna

#mkdir -p $wd/step_10_getFinalGenomeCoordsAndSeqs


#finalBed=$wd/step_10_getFinalGenomeCoordsAndSeqs/${spp}.FinalFunctionalORs.coords.bed
finalBed=$wd/step_10_getFinalGenomeCoordsAndSeqs/finalORCoordinates.all.0based.bed
# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 200 bp on either side (skipping all the stuff Gang did in Step 2 that made sure
# the query start and end matched the hit, instead just giving a big 200 bp buffer. Hopefully works
cd $wd
bedtools getfasta -s -name -bed $finalBed -fi $REFERENCE -fo step_10_getFinalGenomeCoordsAndSeqs/${spp}.step_10_results.finalORntSeqs.all.fna
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 200bp to each side of your range 
