#! /bin/bash
#$ -cwd
#$ -l h_cpu=10:00:00,h_data=10G,highp
#$ -o ./log/otters_20170922/
#$ -e ./log/otters_20170922/
#$ -m abe
#$ -M ab08028
#$ -t 1-300
#$ -N shuffle
# make script executable (only need to run once)
source /u/local/Modules/default/init/modules.sh
module load R/3.4.0
i=$SGE_TASK_ID 
branch=$1 # elut pbra or otters or otters.max50
label=$2 #label for the run; generally paml.dateYouRanPaml.Info  
projectname="${branch}_${label}"

./script/polysel_shufflesets_otters.sh $i $projectname $branch $label

sleep 5m
