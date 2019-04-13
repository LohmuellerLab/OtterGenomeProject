#! /bin/bash
#$ -l h_rt=01:00:00,h_data=500M
#$ -cwd
#$ -N prune_trees
#$ -m abe
#$ -M ab08028
#$ -o guidanceReports
#$ -e guidanceReports
#$ -t 1-77
######## STEP 5a, after GUIDANCE residue masking, removing spp.

# only want to make trees if there are 6 spp remaining after guidance removal of bad spp
# and only want to do it if the species of interest are there (sea otter and giant otter)

source /u/local/Modules/default/init/modules.sh
module load R
res=0.93 # residue masking filter I used 
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file 
allFasta=all15spp_CDS.TranscriptIDs.1.20170606.fna
date=20170608 # date you ran guidance
j=`printf "%02d" $((10#$SGE_TASK_ID - 10#1))`
#j=`printf "%02d" $(($SGE - 1))`
#mkdir group_${j}
start=$((10#$j*200))
echo "start: " $start
end=$((10#$j*200+199))
echo "end: " $end
for (( i=$start; i<=$end; i++ ))
do 
echo $i
p=$(($i+1))
# get human T ID to be overall gene name
TName=`awk -v a=${row1} 'FNR==a {print $"'$p'"}' $clusters`
# but if TName is missing in human, sub in the dog:
if [[ $TName = "NA" ]]
then 
TName=`awk -v a=${row2} 'FNR==a {print $"'$p'"}' $clusters`
fi

cluster=`printf "%05d" $i`
############## need to make trees based on the taxa that are there. #########
rscript=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/scripts/Step5a.2.prune_tree.R
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_20170606 
header=group_${j}/cluster_${cluster}_group_${j}_${TName}
guidanceRun=guidance_${date}
fastaFile=$wd/$header/$guidanceRun/MSA.PRANK.aln.With_Names.ResMask.${res}
fullTree=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_20170606/unrooted_tree_15spp_pandaInRightPlace_horseNotWithCow.USETHIS.20170606.txt

 
# get taxa list for each cluster
grep "^>" $fastaFile | sed 's/>//g' > $wd/$header/$guidanceRun/taxaList.txt
# prune trees:
# usage: Step5a.2.prune_tree.R <fullTree> <header>
# output is called prunedTree.noFG.txt
# NEED TO USE THE TAXALIST IN THE GUIDANCE DIR, NOT IN THE MAIN CLUSTER DIRECTORY:
Rscript $rscript $fullTree $header/$guidanceRun

# then you want to check if elut, pbra, or both are present and make foreground branches:
numTaxa=`wc -l $wd/$header/$guidanceRun/taxaList.txt | awk '{print $1}'`
# must be at least 6 taxa first: 
# don't mind if hsap, cfam not in final alignment, just had to have an ortholog at the beginning
# e.g. human may have more exons, disrupt stuff.
if [ $numTaxa -ge 6 ]
then
if grep -q elut $header/$guidanceRun/prunedTree.noFG.txt
then
sed 's/elut/elut#1/g' $header/$guidanceRun/prunedTree.noFG.txt > $header/prunedTree.elut.post.$guidanceRun.txt
fi

if grep -q pbra $header/$guidanceRun/prunedTree.noFG.txt
then
sed 's/pbra/pbra#1/g' $header/$guidanceRun/prunedTree.noFG.txt > $header/prunedTree.pbra.post.$guidanceRun.txt
fi

if grep -q "elut,pbra" $header/$guidanceRun/prunedTree.noFG.txt
then
sed 's/\(elut,pbra\)/\(elut,pbra\)#1/g' $header/$guidanceRun/prunedTree.noFG.txt > $header/prunedTree.otters.post.$guidanceRun.txt
fi


## optional: add pinnipeds: (only to trees that have BOTH pinnipeds)
#if grep -q "oros,lwed" $header/$guidanceRun/prunedTree.noFG.txt
#then
#sed 's/\(oros,lwed\)/\(oros,lwed\)#1/g' $header/$guidanceRun/prunedTree.noFG.txt > $header/prunedTree.pinnipeds.post.$guidanceRun.txt
#fi


#if grep -q oros $header/$guidanceRun/prunedTree.noFG.txt
#then
#sed 's/oros/oros#1/g' $header/$guidanceRun/prunedTree.noFG.txt > $header/prunedTree.oros.post.$guidanceRun.txt
#fi

#if grep -q lwed $header/$guidanceRun/prunedTree.noFG.txt
#then
#sed 's/lwed/lwed#1/g' $header/$guidanceRun/prunedTree.noFG.txt > $header/prunedTree.lwed.post.$guidanceRun.txt
#fi

#fi

done
sleep 5m
