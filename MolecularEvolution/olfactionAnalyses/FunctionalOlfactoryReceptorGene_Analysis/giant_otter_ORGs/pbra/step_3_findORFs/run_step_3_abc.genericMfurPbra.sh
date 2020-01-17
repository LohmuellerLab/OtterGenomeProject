### run step 3: find ORFs:
spp=$1 # have to use source not sh to run put in mfur or pbra
wd=/work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_3_findORFs
cd $wd
# fasta must be in same directory:
cp /work2/abeichman/gangLi_OlfactionPipeline_${spp}/step_2_getSequencesFromGenome/step_2_results.olfacUniqueBlastHits.stranded.fasta $wd

# have to manually set the file name if it's changed from step2
# a: 
perl step_3_a_Find_ORF-ABmod.pl  # gives you fasta as output

# b: 
perl step_3_b_Find_ORF-ABmod_withPrintStatements.pl # gives you full INFO on each ORF


outputfasta="ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta"
outputinfo="ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta.FULLINFO"
# c: perform some checks
# i. make sure all are > 200 AAs long (600bp)
echo "The following are too short" 
grep "Length:" $outputinfo | awk '{if($8 < 6000)print}' # if anything prints out go and get rid of it

# ii. make sure there are no weird AAs
if grep -q "bad" $outputinfo
then
echo "There are bad codons present"
fi