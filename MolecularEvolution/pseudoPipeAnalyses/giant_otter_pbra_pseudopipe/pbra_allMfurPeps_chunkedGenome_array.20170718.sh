#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=10G
#$ -N pbraCHUNK
#$ -o ppipe_reports_pbra
#$ -e ppipe_reports_pbra
#$ -m abe
#$ -M ab08028
#$ -t 1-100
### PBRA  ##### 
## each group contains a lot of scaffolds ; evenly breaking up genome. See if this works!
source /u/local/Modules/default/init/modules.sh
module load python/2.7
rundate=`date +%Y%m%d`
# need to count from j = 000 to 099

j=`printf "%03d" $((10#$SGE_TASK_ID - 10#1))`

### NEED TO MAKE SURE IT'S ALWAYS A NEW OUTPUT DIR (not sure why...)
# ./pseudopipe.sh [output dir] [masked dna dir] [input dna dir] [input pep dir] [exon dir] 0
ppipe=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes/pseudopipe/bin/pseudopipe.sh
header=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes
mkdir $header/ppipe_output/pbra_allMfurPep_${rundate} # make a overall output dir 

$ppipe \
$header/ppipe_output/pbra_allMfurPep_${rundate}/pteronura_brasiliensis_allMfurPep_out_${j}_${rundate} \
$header/ppipe_input/pteronura_brasiliensis/dna/group_${j}/PteBra.a.lines.15kbOnly.REPEATMASKED.NNs_${j}.fasta \
$header/ppipe_input/pteronura_brasiliensis/dna/group_${j}/%s.fasta \
$header/ppipe_input/pteronura_brasiliensis/allMfurpeps/Mustela_putorius_furo.MusPutFur1.0.pep.all.fa \
$header/ppipe_input/pteronura_brasiliensis/mysql/group_${j}/chr%s_exLocs \
0

sleep 10m
