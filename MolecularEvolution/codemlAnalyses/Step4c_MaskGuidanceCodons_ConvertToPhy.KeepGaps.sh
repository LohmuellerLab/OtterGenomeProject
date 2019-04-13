#! /bin/bash
#$ -l h_rt=10:00:00,h_data=1G,highp,arch=intel*
#$ -cwd
#$ -N maskGuidanceResidues
#$ -m abe
#$ -M ab08028
#$ -o guidanceReports
#$ -e guidanceReports
#$ -t 1-77

source /u/local/Modules/default/init/modules.sh
module load perl
module load ruby

# general params:
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_20170606 
date=20170608
guidanceDir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/guidance.v2.02
guidanceRun=guidance_${date}
trimal=/u/home/a/ab08028/klohmueldata/annabel_data/bin/trimal/source/trimal
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file


res=0.93 # residue masking threshold
gt=0.8 # (gap threshold -- 1 - allowed gap %. so if it is 0.8, then any column with >20% gaps is removed)
## modifying script to replace with lower case x instead of X so you can distinguish downstream.

#### set up 
j=`printf "%02d" $((10#$SGE_TASK_ID - 10#1))`
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
######### Running Guidance:

header="group_${j}/cluster_${cluster}_group_${j}_${TName}"

### STEPS
# 0. write a summary:
# check it hasn't already been masked
if [ ! -e $wd/$header/$guidanceRun/MSA.PRANK.aln.With_Names.ResMask.${res}.forPaml.phy ]
then
echo $date >> $wd/$header/filteringREADME.${date}
echo "Guidance residue masking treshold: $res (no col removed)" >> $wd/$header/filteringREADME.${date}

# Mask AA Residues < res with Guidance
# I modified the Guidance code to replace residues with lower case x and n (to distinguish from any previous X or N, as in Selectome)
# 
perl $guidanceDir/www/Guidance/maskLowScoreResidues.ABedits.x.pl $wd/$header/$guidanceRun/MSA.PRANK.aln.Sorted.With_Names $wd/$header/$guidanceRun/MSA.PRANK.Guidance2_res_pair_res.scr $wd/$header/$guidanceRun/MSA.PRANK.aln.With_Names.ResMask.${res} $res nuc

# Convert masked fasta to paml phylip format
$trimal -phylip_paml -in $wd/$header/$guidanceRun/MSA.PRANK.aln.With_Names.ResMask.${res} -out $wd/$header/$guidanceRun/MSA.PRANK.aln.With_Names.ResMask.${res}.forPaml.phy


fi
done


sleep 5m
