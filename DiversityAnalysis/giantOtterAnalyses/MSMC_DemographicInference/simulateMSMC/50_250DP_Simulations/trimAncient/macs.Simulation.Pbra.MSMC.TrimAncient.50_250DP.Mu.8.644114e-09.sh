#!/bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=28G,highp
#$ -N pbraMSMC1
#$ -m abe
#$ -M ab08028
rundate=`date +%Y%m%d`
model=pbra.msmc.trimAncient.8.644114e-09
mkdir -p ${model}
for j in {1..6}
do
mkdir -p ${model}/group_$j.${model}
cd ${model}/group_$j.${model}
cp ../../../macs ./
cp ../../../msformatter ./
for i in {1..10}
do
# giant otter msmc model -- full msmc
mu=8.644114e-09
r=1e-08
Na=19211.3974344
rho=0.000768455897376
theta=0.000664262038089
date=`date +%Y%m%d`
SEED=$((date+$RANDOM+((j-1)*10)+i))
# this is a new addition! need to have a different random seed for each simulation; if they start within a second of each other, they will have the same seed. not an issue for big simulations of 30Mb because those are slow, but 100kb can start within a second of each other!
./macs 2 30000000 -t 0.000664262038089 -r 0.000768455897376 -s $SEED -eN 0.0 0.108152203196 -eN 0.0110939502447 0.122530155785 -eN 0.0224760699 0.349389439384 -eN 0.034161819732 0.570349897802 -eN 0.046167623982 0.752000599431 -eN 0.058511848896 0.887994266552 -eN 0.071213763978 0.97807274003 -eN 0.084294896877 1.03028374328 -eN 0.097778732301 1.05780054386 -eN 0.111690561474 1.07165204392 -eN 0.12605868648 1.07864709655 -eN 0.140913968634 1.07864709655 -eN 0.15629073174 1.07369659796 -eN 0.17222570829 1.07369659796 -eN 0.18876436227 1.05523136612 -eN 0.20595035115 1.05523136612 -eN 0.22383787041 1.02883678972 -eN 0.24248713725 1.02883678972 -eN 0.26196439059 1.0 -eN 0.28234941822 1.0  > group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt 
./msformatter < group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt
done
cd ../../
done
