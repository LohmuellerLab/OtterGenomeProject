## RUN IN THE SHELL : Step7b.sh elut 20170301 (date you ran paml)
# (this script was formerly named Step6b)
display_usage() { 
	echo "This script requires you to specify species (otters, elut, pbra) and date you ran paml in format YYYYMMDD *and* what replicate it is (1,2,or 3) \n Note that date you aligned orthologous clusters is defined below -- if re-aligned, change it." 
	echo -e "\nUsage:./script.sh [species] [paml run date] [replicate]\n" 
	}
	
# if fewer than three arguments supplied, display usage 
	if [  $# -le 1 ] 
	then 
		display_usage
		exit 1
	fi 
 
# check whether user had supplied -h or --help . If yes display usage 
	if [[ ( $# == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 
date=20170606 # this is the date the clusters were made, not the date you're running this script, 
# eventually should be run by user 
# and this is the date of the codeml run for the output:

flag=$1 # species that is the foreground branch (otters, elut or pbra)
rundate=$2 # date you ran paml; 20170327 for no-swamp; 20170328 for w/ swamp
replicate=$3 # change if a different codeml replicate
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam

echo "species: " $1
echo "paml run on: " $2
echo "replicate: " $3
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_${date}
clusters=$wd/one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt


# dir to gather output files:
mkdir -p $wd/output_codeml_${rundate}
gatherOut=$wd/output_codeml_${rundate}
## set up the header of the output file: 
echo "branch replicate group	cluster	TName ${flag}.lnL_null ${flag}.lnL_alt ${flag}.bebCount sppCount ${flag}.treeLength seqLength_mayInclGaps" > $gatherOut/codeml_${flag}.rep.${replicate}.${rundate}.fullInfo.txt


# this gets you a number that is one less than SGE (because my numbering starts at 0)
for j in {00..77}
do
# start and end for each group of 200 clusters
start=$((10#$j*200))
echo "start: " $start
end=$((10#$j*200+199))
echo "end: " $end
for (( i=$start; i<=$end; i++ ))
do 
#echo $i 
p=$(($i+1))
#echo  "p:" $p
# get human T ID to be overall gene name ; note this is actually cat ID due to addition of pbra 
TName=`awk -v a=${row1} 'FNR==a {print $"'$p'"}' $clusters`
# but if TName is missing in human, sub in the dog:
if [[ $TName = "NA" ]]
then 
TName=`awk -v a=${row2} 'FNR==a {print $"'$p'"}' $clusters`
fi
#echo $TName
cluster=`printf "%05d" $i`
header="group_${j}/cluster_${cluster}_group_${j}_${TName}"
# if didn't use swamp:
#pamlDir=codeml_mafft_noSwamp_model2_PosSelection_cluster_${cluster}_group_${j}_${rundate}
# if used swamp:
pamlDir=codeml_branchSite_results_rep_${replicate}_${rundate}

infileAlt=codemlResult_${cluster}_group_${j}_${TName}_${rundate}.altModel.$flag.rep.${replicate}.mlc

infileNull=codemlResult_${cluster}_group_${j}_${TName}_${rundate}.nullModel.$flag.rep.${replicate}.mlc
# first check if there is a result: 
if [ -e $wd/$header/$pamlDir/$infileAlt ]
then
#echo $sppCount
treeLength=`grep "tree length" $wd/$header/$pamlDir/$infileNull | awk '{print $4}'`
#echo $treeLength
# need to do things differently if gaps were removed or not (mlc has diff structure)
if grep --quiet "gaps" $wd/$header/$pamlDir/$infileNull
then
# if gaps:
len=""
sppCount=""
len=`grep "After deleting gaps." $wd/$header/$pamlDir/$infileNull | awk '{print $4}'`
sppCount=`head -n3 $wd/$header/$pamlDir/$infileNull | tail -n1 | awk '{print $1}'`

else
# if no gaps:
len=""
sppCount=""
len=`head -n1 $wd/$header/$pamlDir/$infileNull | awk '{print $2}'`
sppCount=`head -n1 $wd/$header/$pamlDir/$infileNull | awk '{print $1}'`

fi

lnL_null=""
lnL_alt=""
lnL_null=`grep lnL $wd/$header/$pamlDir/$infileNull | awk '{print $5}'`
lnL_alt=`grep lnL $wd/$header/$pamlDir/$infileAlt | awk '{print $5}'`

# BEB sites
mkdir -p $wd/$header/$pamlDir/bebInfo
# want to get all BEB codons, even non-significant.
grep -h "Bayes Empirical Bayes (BEB) analysis" -A 5000 $wd/$header/$pamlDir/$infileAlt | grep "The grid" -B 5000 | grep -E "[0-9]+ [A-Z]+" > $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEB.sites.txt

# get selected codons:
codons=""
codons=`awk -vORS="\t" '{print $1}' $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEB.sites.txt | sed 's/ $/\n/g'`
# empty file:
> $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEBcodons.txt
# get total count of BEBs that are significant  
bebCount=""
bebCount=`grep -c "\*" $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEB.sites.txt`
# get codons under selection: 
# actually want to look at all BEB sites?
for k in $codons
do
if grep --quiet "gaps" $wd/$header/$pamlDir/$infileAlt
then
echo "there are gaps"
# gaps:
grep -E "^(\s)+$k" $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEB.sites.txt >> $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEBcodons.txt
grep "After deleting gaps." $wd/$header/$pamlDir/$infileAlt -A $((sppCount+2)) | tail -n$sppCount | awk '{FS=" "; print $1,$"'$(($k+1))'"}' >> $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEBcodons.txt
else 
# no gaps:
echo "there aren't gaps"
grep -E "^(\s)+$k" $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEB.sites.txt >> $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEBcodons.txt
head -n $((sppCount+2)) $wd/$header/$pamlDir/$infileAlt | head -n$((sppCount+2)) | tail -n$sppCount | awk '{FS=" "; print $1,$"'$(($k+1))'"}' >> $wd/$header/$pamlDir/bebInfo/codemlResult_${cluster}_group_${j}_${TName}_${rundate}.${flag}.rep.${replicate}.sigBEBcodons.txt
fi
done

## check if any variable is empty and set to NA
if [ -z "$lnL_null" ]
then
lnL_null="NA"
fi
if [ -z "$lnL_alt" ]
then
lnL_alt="NA"
fi
if [ -z "$bebCount" ]
then
bebCount="NA"
fi
if [ -z "$sppCount" ]
then
sppCount="NA"
fi
if [ -z "$treeLength" ]
then
treeLength="NA"
fi
if [ -z "$len" ]
then
len="NA"
fi
# then echo them all: 
echo "$flag	$replicate	$j	$cluster	$TName	$lnL_null	$lnL_alt	$bebCount	$sppCount	$treeLength	$len" >> $gatherOut/codeml_${flag}.rep.${replicate}.${rundate}.fullInfo.txt

#fi
else 
echo "there is no $infileAlt in $header"
fi
cd $wd
done
done


