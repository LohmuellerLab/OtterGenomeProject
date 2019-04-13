#! /bin/bash
#$ -l h_rt=100:00:00,h_data=2G,highp
#$ -cwd
#$ -pe shared 5
#$ -N guidance
#$ -m abe
#$ -M ab08028
#$ -o guidanceReports
#$ -e guidanceReports
#$ -t 1-77
# usage: qsub Step4-alt.sh <todays date 20170608>


source /u/local/Modules/default/init/modules.sh
module load perl
module load ruby

row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_20170606
date=$1 # date you're running guidance
boot=30 #number of bootstraps; prank is slow so only doing 30 (in Guidance manual)
seqCutoff=0.6 # default:0.6 sequence cutoff to drop bad sequences (species) guidance default.
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
# but if TName is missing in human, sub in the dog name:
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

# Guidance USAGE: perl $guidanceDir/www/Guidance/guidance.pl --seqFile SEQFILE --msaProgram [MAFFT|PRANK|CLUSTALW|MUSCLE] --seqType [aa|nuc|codon] --outDir FULL_PATH_OUTDIR
# --MSA_Param: passing parameters for the alignment program e.g -F to prank. To pass parameter containning '-' in it, add \ before each '-' e.g. \-F for PRANK

#  --proc_num: number of processors to use (default=1)
# --seqCutoff: Confidence cutoff between 0 to 1. Default=0.6
# --colCutoff: Confidence cutoff between 0 to 1. Default=0.93
# --outOrder [aligned|as_input] default=aligned
 
# --seqFile: Input sequence file in FASTA format
#--msaProgram: Which MSA program to use
#--seqType: Type of sequences for alignment (amino acids, nucleotides, or codons)
#--outDir: Output directory that will be created automatically and hold all output files [please provid full (and not relative) path]

# Prank codon +F:
## need full path to sequence and to output
perl $guidanceDir/www/Guidance/guidance.pl \
--seqFile $wd/$header/$fasta \
--msaProgram PRANK --prank $prank \
--seqType codon \
--outDir $wd/$header/guidance_${date} \
--bootstraps $boot --seqCutoff $seqCutoff --colCutoff $colCutoff \
--proc_num 5 \
--outOrder as_input \
--MSA_Param '\-F'
# MSA_Param F is only for prank. 

# the outputs that matter are:
# MSA.PRANK.aln.With_Names this is the alignment without bad columns removed
# MSA.PRANK.Without_low_SP_Col.With_Names this is the alignment with bad columns removed
# Seqs.Orig_DNA.fas.FIXED.Without_low_SP_Seq.With_Names this is the original unaligned seqs
# with bad species removed. You would re-run Guidance on this. 
# how to check if species were removed? 
# MSA.PRANK.Guidance2_res_pair_seq.scr_with_Names shows the scores of the sequences (could check if any are < 0.6)
# 1. see if any species were removed (below 0.6)
# check how frequently this happens in a group. 
# 2. make sure that pbra or elut is still present
# 3. rerun guidance once bad sequences are removed

done
sleep 5m

