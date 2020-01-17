#### run on SIRIUS:

# deduped ref genome:
spp="mfur"
wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}


REFERENCE=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_1_ab_blastDb_blastResults/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 

finalBed=$wd/step_10_getFinalGenomeCoordsAndSeqs/finalORCoordinates.all.0based.bed
# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 300 bp on either side
cd $wd
bedtools getfasta -s -name -bed $finalBed -fi $REFERENCE -fo step_10_getFinalGenomeCoordsAndSeqs/${spp}.step_10_results.finalORntSeqs.all.fna
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 300bp to each side of your range 
