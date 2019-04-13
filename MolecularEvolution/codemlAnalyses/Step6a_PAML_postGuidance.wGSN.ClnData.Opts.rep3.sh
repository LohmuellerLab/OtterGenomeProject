#! /bin/bash
#$ -l h_rt=150:00:00,h_data=16G,highp
#$ -cwd
#$ -N paml3
#$ -m abe
#$ -M ab08028
#$ -o pamlReports
#$ -e pamlReports
#$ -t 1-77
# run the third replicate of codeml
display_usage() {
	echo "This is for rep 3 (omega=2 starting). This script requires you to specify: 1) species for the foreground branch (otters, elut, pbra) \n 2) rundate for this run of paml (YYYYMMDD) \n 3) G/S/N (gblocks, swamp, none) 4) (optional) if S, need to specify swamp date \n Note that date you aligned orthologous clusters is defined below -- if re-aligned, change it." 
	echo -e "\nUsage: qsub script.sh [species] [date] [filter: G/GMED/S/N (gblocks, swamp or none)] [CLEANDATA 0 1] OPTIONAL if S: swamp run date \n" 
	}
	
# if fewer than three arguments supplied, display usage 
	if [  $# -lt 4 ] 
	then 
		display_usage
		exit 1
	fi 
# if swamp was run, but no swamp date given:
	if [ $3 = "S" ]
	then
	if [ $5 = "" ]
	then
	echo "You need to specify swamp rundate as fourth input"
	exit 1
	fi
	fi 
# check whether user had supplied -h or --help . If yes display usage 
	if [[ ( $# == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 

# change this to be whether you want data cleaned or not:
clnData=$4 # 1 for clean, 0 for don't clean
### to submit: qsub Step6_generic.blahblah.sh elut (or otters or pbra)
res=0.93 # this is whatever Guidance residue masking filter you used in the previous step
rep=3 # if you're running paml multiple times for convergence  
# allows you to try different guidance filters 
date=20170606 # the date your clusters were made, not the date you run this script
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_${date}
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file 
# and this is the date of the codeml run for the output:
#rundate=`date +%Y%m%d` # if you think the script will start after the next day, set by hand
rundate=$2 # set it in the script
# guidanceRun: this is the guidance run you're using
guidanceRun=guidance_20170608 # change this if you redo Guidance
flag=$1 # this is the branch (elut, pbra, otters, pinnipeds, oros, lwed)
codeml=/u/home/a/ab08028/klohmueldata/annabel_data/bin/paml4.9c/bin/codeml
# generic ctl file for null model (omega fixed to 1) CLEANDATA=0 (won't remove ambiguities)
genericNull=$wd/paml_generic_files/generic.codeml.fixed.null.CleanData${clnData}.ctl
# generic ctl file for alt model (omega allowed to vary)  CLEANDATA=0 (won't remove ambiguities)
genericAlt=$wd/paml_generic_files/generic.codeml.alt.CleanData${clnData}.omega2.rep3.ctl
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
# got the whitespace correct here:
if [[ $TName = "NA" ]] 
then 
TName=`awk -v a=${row2} 'FNR==a {print $"'$p'"}' $clusters`
fi
cluster=`printf "%05d" $i`
header="group_${j}/cluster_${cluster}_group_${j}_${TName}"

####### have to check if correct tree is present, else move on.
if [ -e $wd/$header/prunedTree.${flag}.post.${guidanceRun}.txt ]
then
echo "$flag tree is present"

tree=$wd/$header/prunedTree.${flag}.post.${guidanceRun}.txt
# note that the version on HOffman doesn't have the "rep" label for alt models -- fix in future runs
outfileNull=codemlResult_${cluster}_group_${j}_${TName}_${rundate}.nullModel.$flag.rep.$rep.mlc
outfileAlt=codemlResult_${cluster}_group_${j}_${TName}_${rundate}.altModel.$flag.rep.$rep.mlc
# go to directory: 
cd $wd/$header
# make a directory for this paml run:
mkdir -p codeml_branchSite_results_rep_${rep}_${rundate}
pamlDir=codeml_branchSite_results_rep_${rep}_${rundate}

# check if paml has already been run with these settings
# check alt file because it only is made once null model is finished
if ! grep -q "lnL" $pamlDir/$outfileAlt
then
# check if gblocks or swamp was used
filter=$3 # options are G (gblocks), S (Swamp), N (none)
# single = is ok. 
### WHITE SPACE MATTERS@!!!! [ $filter="G" ] WILL ALWAYS EVALUATE AS TRUE. [ $filter = "G" ] IS CORRECT
#(that's insane.)

if [ $filter = "G" ]
then
echo "filter was Gblocks"
seqName=MSA.PRANK.aln.With_Names.ResMask.${res}-gb.phylip # this is w/ gblocks! 
cp $wd/$header/$guidanceRun/$seqName $pamlDir

elif [ $filter = "GMED" ]
then
echo "filter was Gblocks-medium (block size 30)"
seqName=MSA.PRANK.aln.With_Names.ResMask.${res}gbmed.phylip # this is w/ gblocks medium
cp $wd/$header/$guidanceRun/$seqName $pamlDir

elif [ $filter = "N" ]
then
echo "filter was None"
seqName=MSA.PRANK.aln.With_Names.ResMask.${res}.forPaml.phy # this is without gblocks or swamp
cp $wd/$header/$guidanceRun/$seqName $pamlDir
elif [ $filter = "S" ]
then
swampDate=$5 ### change if you redo swamp!
if [ $swampDate = "" ]
then
echo "Need to specify swamp date"
exit 1
fi
echo "filter was Swamp $swampDate"
codemlDir=swamp_codeml0_${swampDate}
seqName=MSA.PRANK.aln.With_Names.ResMask.${res}.forPaml.mask1_and_mask2.phy
cp $wd/$header/$codemlDir/$seqName $pamlDir
else
echo "FILTER MUST BE G, GMED, S or N(upper case)"
fi

# need to put ctl file and ## generate null control file :
cp $genericNull $pamlDir/codeml_fixed.null.${flag}.ctl
echo "seqfile = " $seqName >> $pamlDir/codeml_fixed.null.${flag}.ctl
echo "treefile = " $tree >> $pamlDir/codeml_fixed.null.${flag}.ctl
echo "outfile = " $outfileNull >> $pamlDir/codeml_fixed.null.${flag}.ctl
 
## generate alternative control file 
cp $genericAlt $pamlDir/codeml_alt.${flag}.ctl
echo "seqfile = " $seqName >> $pamlDir/codeml_alt.${flag}.ctl
echo "treefile = " $tree >> $pamlDir/codeml_alt.${flag}.ctl
echo "outfile = " $outfileAlt >> $pamlDir/codeml_alt.${flag}.ctl

# get masked .phy file in there:

# go to paml dir and run codeml:
cd $pamlDir
# run null
$codeml ./codeml_fixed.null.${flag}.ctl
# run alt
$codeml ./codeml_alt.${flag}.ctl

fi
cd $wd 
fi

done

sleep 5m

