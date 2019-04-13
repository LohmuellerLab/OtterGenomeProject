#! /bin/bash
#$ -cwd
#$ -l h_cpu=5:00:00,h_data=6G,highp
#$ -o ./log/otters_20170922/
#$ -e ./log/otters_20170922/
#$ -m abe
#$ -M ab08028
#$ -N allTests
# make script executable (only need to do once)
# usage script <branch> <date>
source /u/local/Modules/default/init/modules.sh
module load R/3.4.0
branch=$1 # elut pbra or otters
label=$2 #label of paml run; date you ran paml 
projectname="${branch}_${label}"

./script/polysel_cluster.otters.sh $projectname $branch $label
