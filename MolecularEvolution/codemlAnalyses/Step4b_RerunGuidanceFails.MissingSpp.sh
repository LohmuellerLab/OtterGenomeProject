#! /bin/bash
#$ -l h_rt=50:00:00,h_data=4G,highp,arch=intel*
#$ -cwd
#$ -pe shared 5
#$ -N guidanceRedos
#$ -m abe
#$ -M ab08028
#$ -o guidanceReports
#$ -e guidanceReports
#$ -t 1-77
## Check a series of errors and re-run guidance if any alignments failed:

source /u/local/Modules/default/init/modules.sh
module load perl
module load ruby
seqtk=/u/home/a/ab08028/klohmueldata/annabel_data/bin/seqtk/seqtk
iteration=2 # leave this empty if it's the first run of RedoFails, or as 2 if you're running again (to catch the last failures)
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_20170606
date=20170608 ## date you're running guidance; be careful of using the same dir, will break Guidance
boot=30 #number of bootstraps; prank is slow so only doing 30 (in Guidance manual)
seqCutoff=0.6 # default:0.6 sequence cutoff to drop bad sequences (species)
colCutoff=0.93 # default: 0.93 (paper that tested Guidance went up to 0.99) drop columns below this consistency 
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
guidanceDir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/guidance.v2.02
header="group_${j}/cluster_${cluster}_group_${j}_${TName}"
fasta="Seqs_cluster_${cluster}_group_${j}_${TName}.fna"
prank=/u/home/a/ab08028/klohmueldata/annabel_data/bin/prank/bin/prank


## Check if there are internal stop codons:
# see if the guidance log shows stop codon problem
if grep -q "A Stop" $wd/$header/guidance_${date}/log
then
badSppNum=`grep "A Stop" $wd/$header/guidance_${date}/log | awk '{print $9}'`
badSppID=`grep $badSppNum $wd/$header/guidance_${date}/Seqs.Codes | awk '{print $1}'`
# if don't have one already, make a preguidance taxa list:
if [ ! -e $wd/$header/taxaList.preGuidance.txt ]
then
grep "^>" $wd/$header/$fasta | sed 's/>//g' > $wd/$header/taxaList.preGuidance.txt
fi
# move old fasta:
cp $wd/$header/$fasta $wd/$header/$fasta.CONTAINSBADSPP.${iteration}
# make subseq taxa list:
grep -v $badSppID $wd/$header/taxaList.preGuidance.txt > $wd/$header/taxaList.RemovedBadSpp.${iteration}.txt
$seqtk subseq $wd/$header/$fasta.CONTAINSBADSPP.${iteration} $wd/$header/taxaList.RemovedBadSpp.${iteration}.txt > $wd/$header/$fasta
# re run guidance:
echo "There was a sequence with a premature stop codon" >> $wd/$header/Guidance.${date}.${iteration}.log
echo $header >> $wd/GuidanceSecondRound.StopCodon.${date}.${iteration}.log
mv $wd/$header/guidance_${date} $wd/$header/guidance_${date}_orig${iteration}_FAILED
perl $guidanceDir/www/Guidance/guidance.pl \
--seqFile $wd/$header/$fasta \
--msaProgram PRANK --prank $prank \
--seqType codon \
--outDir $wd/$header/guidance_${date} \
--bootstraps $boot --seqCutoff $seqCutoff --colCutoff $colCutoff \
--proc_num 5 \
--outOrder as_input \
--MSA_Param '\-F'
fi

## check if sequences have a div by 3 problem:

if grep -q "not divisible by 3" $wd/$header/guidance_${date}/log
then
badSppNum=`grep "not divisible by 3" $wd/$header/guidance_${date}/log | awk '{print $3}'`
badSppID=`grep $badSppNum $wd/$header/guidance_${date}/Seqs.Codes | awk '{print $1}'`
# if don't have one already, make a preguidance taxa list:
if [ ! -e $wd/$header/taxaList.preGuidance.txt ]
then
grep "^>" $wd/$header/$fasta | sed 's/>//g' > $wd/$header/taxaList.preGuidance.txt
fi
# move old fasta:
cp $wd/$header/$fasta $wd/$header/$fasta.CONTAINSBADSPP.${iteration}
# make subseq taxa list:
grep -v $badSppID $wd/$header/taxaList.preGuidance.txt > $wd/$header/taxaList.RemovedBadSpp.${iteration}.txt
$seqtk subseq $wd/$header/$fasta.CONTAINSBADSPP.${iteration} $wd/$header/taxaList.RemovedBadSpp.${iteration}.txt > $wd/$header/$fasta
# re run guidance:
echo "There was a sequence not divisible by 3" >> $wd/$header/Guidance.${date}.${iteration}.log
echo $header >> $wd/GuidanceSecondRound.DivBy3.${date}.${iteration}.log
mv $wd/$header/guidance_${date} $wd/$header/guidance_${date}_orig${iteration}_FAILED
perl $guidanceDir/www/Guidance/guidance.pl \
--seqFile $wd/$header/$fasta \
--msaProgram PRANK --prank $prank \
--seqType codon \
--outDir $wd/$header/guidance_${date} \
--bootstraps $boot --seqCutoff $seqCutoff --colCutoff $colCutoff \
--proc_num 5 \
--outOrder as_input \
--MSA_Param '\-F'
fi


### check If guidance failed randomly (generally a memory error)
# if there was an error directory, try rerunning Guidance:
if [ -d $wd/$header/guidance_${date}_err ]
then
echo "Guidance failed the first time. Trying again" >> $wd/$header/Guidance.${date}.${iteration}.log
echo $header >> $wd/GuidanceSecondRound.RedoFails.${date}.${iteration}.log
mv $wd/$header/guidance_${date} $wd/$header/guidance_${date}_orig${iteration}_FAILED
mv $wd/$header/guidance_${date}_err $wd/$header/guidance_${date}_err_orig${iteration}_FAILED
perl $guidanceDir/www/Guidance/guidance.pl \
--seqFile $wd/$header/$fasta \
--msaProgram PRANK --prank $prank \
--seqType codon \
--outDir $wd/$header/guidance_${date} \
--bootstraps $boot --seqCutoff $seqCutoff --colCutoff $colCutoff \
--proc_num 5 \
--outOrder as_input \
--MSA_Param '\-F'
fi


### check if a species was flagged as disrupting the alignment
# check if the species were filtered, if they were, rerun Guidance on the new input
# if file isn't empty, want to rerun guidance 
# the Removed_Seq file shows what sequences were taken out; if empty then you don't need to do antyhing
# if not empty, then rerun using Seqs.Orig_DNA.fas.FIXED.Without_low_SP_Seq.With_Names as input
if [ -e $wd/$header/guidance_${date}/Seqs.Orig_DNA.fas.FIXED.Removed_Seq.With_Names ]
then
# need to move guidance folder
echo "Guidance flagged some species to be removed. Re-running Guidance without these spp:" >> $wd/$header/Guidance.${date}.${iteration}.log
grep "^>" $wd/$header/guidance_${date}/Seqs.Orig_DNA.fas.FIXED.Removed_Seq.With_Names >> $wd/$header/Guidance.${date}.${iteration}.log
echo $header >> $wd/GuidanceSecondRound.RemoveSpp.${date}.${iteration}.log
mv $wd/$header/guidance_${date} $wd/$header/guidance_${date}_orig${iteration}_allInputSpp
perl $guidanceDir/www/Guidance/guidance.pl \
--seqFile $wd/$header/guidance_${date}_orig${iteration}_allInputSpp/Seqs.Orig_DNA.fas.FIXED.Without_low_SP_Seq.With_Names \
--msaProgram PRANK --prank $prank \
--seqType codon \
--outDir $wd/$header/guidance_${date} \
--bootstraps $boot --seqCutoff $seqCutoff --colCutoff $colCutoff \
--proc_num 5 \
--outOrder as_input \
--MSA_Param '\-F'
fi

done


sleep 5m
