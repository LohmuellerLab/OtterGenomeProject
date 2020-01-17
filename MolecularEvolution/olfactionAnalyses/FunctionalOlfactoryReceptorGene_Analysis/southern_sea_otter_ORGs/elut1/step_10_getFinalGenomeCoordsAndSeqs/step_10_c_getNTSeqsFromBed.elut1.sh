#### run on SIRIUS:

# deduped ref genome:
spp=elut1
wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}

REFERENCE=/data3/abeichman/gidgetDeNovoGenome/DovetailAssembly/sea_otter_23May2016_bS9RH.fasta



finalBed=$wd/step_10_getFinalGenomeCoordsAndSeqs/${spp}.FinalFunctionalORs.coords.bed
# this should be the 0based start and stop coordinates of each unique match (cleaned in step 1c)
# plus and minus 300 bp on either side
cd $wd
bedtools getfasta -s -name -bed $finalBed -fi $REFERENCE -fo step_10_getFinalGenomeCoordsAndSeqs/${spp}.step_10_results.finalORntSeqs.fna
# use -s so that it is stranded! will pull reverse complement as needed
# use the name column to name it
# make sure you've added +- 300bp to each side of your range 
