#! /bin/bash
#$ -l h_rt=10:00:00,h_data=10G,highp
#$ -cwd
#$ -N swampPostGuidance
#$ -m abe
#$ -M ab08028
#$ -o swampReports
#$ -e swampReports
#$ -t 1-77
############ Swamp prep script:
# run this after guidance 
source /u/local/Modules/default/init/modules.sh
module load python 
module load perl
module load R/3.4.0
####*GUIDANCE *####
dirDate=20170606 # this is the date the clusters were made, not the date you're running this script
rundate=20170713 #`date +%Y%m%d` # and this is the date of the codeml run for the output:
GuidDate=20170608 # date you ran guidance
res=0.93 # residue masking filter I used 
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam

### programs:
codeml=/u/home/a/ab08028/klohmueldata/annabel_data/bin/paml4.9c/bin/codeml # path to codeml
swamp=/u/home/a/ab08028/klohmueldata/annabel_data/bin/SWAMP/SWAMP.py
#### paths to scripts/files you'll be using: (swampMaterials dir with branches, scripts etc)
swampDir=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/swampMaterials_20170228
generic0ctl=$swampDir/generic.Model0.forSwamp.clean0.ctl  # path to generic 0 ctl model
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file 
# for cluster names

#### swamp paramaters ###
t1=2 # AA changes within window 
w1=3 # window size in codons
m1=40 # min length 40 codons, 120bp
t2=5 # second swamp round
w2=15 # second swamp round 
onlyOtters="Y" # whether you want elut and pbra branches masked, or all

########## get the loop going ##########
# this gets you a number that is one less than SGE (because my numbering starts at 0)
j=`printf "%02d" $((10#$SGE_TASK_ID - 10#1))`
# start and end for each group of 200 clusters

start=$((10#$j*200))
echo "start: " $start
end=$((10#$j*200+199))
echo "end: " $end
for (( i=$start; i<=$end; i++ ))
do 
#echo $i
p=$(($i+1))
# get human T ID to be overall gene name; note: this actually gets CAT id. it's okay though
TName=`awk -v a=${row1} 'FNR==a {print $"'$p'"}' $clusters`
# but if TName is missing in human, sub in the dog:
if [[ $TName = "NA" ]]
then 
TName=`awk -v a=${row2} 'FNR==a {print $"'$p'"}' $clusters`
fi
cluster=`printf "%05d" $i`
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_${dirDate}
cd $wd
header=group_${j}/cluster_${cluster}_group_${j}_${TName}
guidanceRun=guidance_${GuidDate}

########### file names ########################
# if running BEFORE gblocks: this is if you're converting fasta to phylip (before gblocks)
inputPhyName=MSA.PRANK.aln.With_Names.ResMask.${res}.forPaml.phy
#outputPhyName=${header}/$guidanceRun/MSA.PRANK.aln.With_Names.ResMask.${res}.SWAMP.${rundate}.phy

# if running AFTER gblocks:
#inputFasName=PrankAlignment_CODONMODEL_${date}_cluster_${cluster}_group_${j}_HsapID_${TName}.best.fas-gb.NOSPACES
#outputPhyName=PrankAlignment_CODONMODEL_${date}_cluster_${cluster}_group_${j}_HsapID_${TName}.best.sorted.gblocks.phy


######### go to cluster dir:
clusterDir=$header
cd $clusterDir
# make dir to run codeml in:
mkdir -p swamp_codeml0_${rundate}
codemlDir=swamp_codeml0_${rundate}


# Need to sort mafft alignment file alphabetically, and convert to phylip format
# the phylip format needs to have an return after each spp name (not default in .pl script so do it with sed)

if [ -e $wd/$header/$guidanceRun/$inputPhyName ]
then
if [ ! -e $codemlDir/${inputPhyName%.phy}.mask1_and_mask2.phy ]
then
# move the file in there and put an enter between spp name and sequence | and get rid of emptylines
sed -E 's/^([a-z]+)( )+/\1\n/g' $wd/$header/$guidanceRun/$inputPhyName | sed '/^\s*$/d'  >  $codemlDir/$inputPhyName # 5b. Need to run codeml model 0 on it
# Generate ctl file:
## generate null control file 
cd $codemlDir
outfileModel0=codemlResult_cluster_${cluster}_group_${j}_${TName}_${rundate}.Model0.mlc # name of outfile
specificCtl=codeml_${cluster}_group_${j}_${TName}_${rundate}.Model0.swamp.ctl # this cluster's CTL file
tree=$wd/$header/prunedTree.noFG.txt # this is post-guidance even though it doesn't say so 
cp $generic0ctl $specificCtl
echo "seqfile = " $inputPhyName >> $specificCtl
echo "treefile = " $tree >> $specificCtl
echo "outfile = " $outfileModel0 >> $specificCtl
# to get branch labels:

# run codeml:
$codeml $specificCtl
# after codeml get branch names from rst file
### RUN R SCRIPT TO GET BRANCH LABELS ### 
getBranchLabels=$swampDir/getBranchLabelsForSwamp.R
Rscript $getBranchLabels $wd/$header/$codemlDir # this outputs a file called branchLabels.specific.txt
branch_all=$wd/$header/$codemlDir/branchLabels.specific.txt
# don't get rid of spaces
#sed -i'' 's/, /,/g' $branch_all
# if you only want to mask otter branches: 
if [ $onlyOtters==Y ]
then
## pull out the elut and pbra lines 
awk '{if($2=="elut")print}' $branch_all > $wd/$header/$codemlDir/branchLabels.specific.OttersOnly.txt
awk '{if($2=="pbra")print}' $branch_all >> $wd/$header/$codemlDir/branchLabels.specific.OttersOnly.txt
awk '{if($2=="elut, pbra")print}' $branch_all >> $wd/$header/$codemlDir/branchLabels.specific.OttersOnly.txt
awk '{if($2=="pbra, elut")print}' $branch_all >> $wd/$header/$codemlDir/branchLabels.specific.OttersOnly.txt
branch=$wd/$header/$codemlDir/branchLabels.specific.OttersOnly.txt
else
branch=$branch_all
fi

# then want to go up a level
cd ../
# and run swamp on that folder:
### don't use -s for now! blocks out too much sequence w/out gblocks
# needs phylip file to have enters after each name 
$swamp -i $codemlDir -b $branch -t $t1 -w $w1 -m $m1

echo "Running swamp: $rundate. Input: $inputPhyName \n Parameters: round 1: $t1 / $ w1 ; round 2: $t2 / $w2 ; no S." > $wd/$header/swampParameters.${rundate}.txt
if [ $onlyOtters==Y ]
then
echo "only on otter branches" >> $wd/$header/swampParameters.${rundate}.txt
fi

# then want to run swamp again so change name of output file 
# rename output:
mv $codemlDir/${inputPhyName%.phy}_masked.phy $codemlDir/${inputPhyName%.phy}.mask1.phy
# run swamp again:
$swamp -i $codemlDir -b $branch -t $t2 -w $w2 -m $m1
mv $codemlDir/${inputPhyName%.phy}.mask1_masked.phy $codemlDir/${inputPhyName%.phy}.mask1_and_mask2.phy
mv $codemlDir/${inputPhyName%.phy}_masked.phy $codemlDir/${inputPhyName%.phy}.mask2.phy

fi
fi
cd $wd # go back to wd and start again 
done

sleep 5m
