#! /bin/bash
#$ -cwd
#$ -l h_rt=50:00:00,h_data=10G,highp
#$ -N mfurPpipe
#$ -o ppipe_reports_mfur
#$ -e ppipe_reports_mfur
#$ -m abe
#$ -M ab08028
#$ -t 1-84

#### mfur: 84 chunks
# to set up genome (set it up on sirius)
# set up: /work2/abeichman/prepMFUR_forPpipe on Sirius
source /u/local/Modules/default/init/modules.sh
module load python/2.7
rundate=`date +%Y%m%d`
j=`printf "%03d" $((10#$SGE_TASK_ID - 10#1))`

### NEED TO MAKE SURE IT'S ALWAYS A NEW OUTPUT DIR 
# ./pseudopipe.sh [output dir] [masked dna dir] [input dna dir] [input pep dir] [exon dir] 0
ppipe=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes/pseudopipe/bin/pseudopipe.sh
header=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pseudoPipe/pgenes
mkdir $header/ppipe_output/mfur_chunked_MfurPep_${rundate} # make a overall output dir 
$ppipe \
$header/ppipe_output/mfur_chunked_MfurPep_${rundate}/mfur_out_chunk_${j}_${rundate} \
$header/ppipe_input/mustela_putorius_furo/dna/group_${j}/Mustela_putorius_furo.MusPutFur1.0.dna_rm.toplevel.RepeatMaskedByEnsembl.15kbOnly_${j}.fasta \
$header/ppipe_input/mustela_putorius_furo/dna/group_${j}/%s.fasta \
$header/ppipe_input/mustela_putorius_furo/pep/Mustela_putorius_furo.MusPutFur1.0.pep.all.fa \
$header/ppipe_input/mustela_putorius_furo/mysql/group_${j}/chr%s_exLocs \
0

## use ferret peptides 
sleep 10m


