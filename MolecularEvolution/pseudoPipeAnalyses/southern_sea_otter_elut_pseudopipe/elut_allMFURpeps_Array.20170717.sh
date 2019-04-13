#! /bin/bash
#$ -cwd
#$ -l h_rt=330:00:00,h_data=10G,highp
#$ -N pseudoOlfactElut
#$ -o ppipe_reports_mfur
#$ -e ppipe_reports_mfur
#$ -m abe
#$ -M ab08028
#$ -t 1-10
### ELUT ##### 
### use: prepElutforPseudoPipe_chunked.20170320 to set up chunked genome 
source /u/local/Modules/default/init/modules.sh
module load python/2.7
rundate=`date +%Y%m%d`
j=`printf "%02d" $((10#$SGE_TASK_ID - 10#1))`

### NEED TO MAKE SURE IT'S ALWAYS A NEW OUTPUT DIR 
# ./pseudopipe.sh [output dir] [masked dna dir] [input dna dir] [input pep dir] [exon dir] 0
ppipe=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes/pseudopipe/bin/pseudopipe.sh
header=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes
mkdir $header/ppipe_output/elut_chunked_MfurPep_${rundate} # make a overall output dir 
$ppipe \
$header/ppipe_output/elut_chunked_MfurPep_${rundate}/enhydra_lutris_out_chunk_${j}_${rundate} \
$header/ppipe_input/enhydra_lutris_chunked/dna/group_${j}/sea_otter_23May2016_bS9RH.15kbONLY.RepeatRunner.RepeatMasker.MASKED_${j}.fasta \
$header/ppipe_input/enhydra_lutris_chunked/dna/group_${j}/%s.fasta \
$header/ppipe_input/enhydra_lutris_chunked/allMfurpeps/Mustela_putorius_furo.MusPutFur1.0.pep.all.fa \
$header/ppipe_input/enhydra_lutris_chunked/mysql/group_${j}/chr%s_exLocs \
0
## use ferret peptides 
sleep 10m
