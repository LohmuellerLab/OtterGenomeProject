#! /bin/bash
#$ -l h_rt=00:10:00,h_data=100M,highp
#$ -cwd
#$ -N sort_transcripts
#$ -m abe
#$ -M ab08028
#$ -o out1.txt
#$ -e err1.txt
#$ -t 1-77
# 15318 genes, so that is 77 groups 

row1=11 # row of ortholog cluster table where hsap ID is
row2=8 # alt cluster name if there no hsap ID, use cfam
clusters="one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt" # clusters file 
allFasta=all15spp_CDS.TranscriptIDs.1.20170606.fna # fasta file with all cds sequences with transcript IDs

# set up folder numbers:
# this adds padding and subtracts 1 from task id
# need to add 10# to make sure it interprets at base 10 not octal, this works
# want to have more folders with fewer than 1000. 200 per folder. So an array of 47

j=`printf "%02d" $((10#$SGE_TASK_ID - 10#1))`
mkdir group_${j}
start=$((10#$j*200))
echo "start: " $start
end=$((10#$j*200+199)) 
echo "end: " $end
for (( i=$start; i<=$end; i++ ))
do 
echo $i
p=$(($i+1))
# get human T ID to be overall gene name
# with dog as a backup for those with hsap = NA
TName=`awk -v a=${row1} 'FNR==a {print $"'$p'"}' $clusters`
if [[ $TName = "NA" ]]
then 
TName=`awk -v a=${row2} 'FNR==a {print $"'$p'"}' $clusters`
fi
cluster=`printf "%05d" $i`
mkdir group_${j}/cluster_${cluster}_group_${j}_${TName}
# also want to get leading zeroes; use printf: %0d says print pad so that final output is 5 digits long(?)
# so if i is 1999, $cluster will be 01999 (helps keep files in order)
awk '{print $"'$p'"}' $clusters > group_${j}/cluster_${cluster}_group_${j}_${TName}/TrIDS_cluster_${cluster}_group_${j}_${TName}.txt
done 

## This works, but makes extra folders for the entries above 11248, need to manually
# delete those that end without an ID:
#ls group_11/*HsapID_ 
#rm -r group_11/*HsapID_

#so group 11 only has 248 in it

sleep 5m
