#!/bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=28G,highp
#$ -N elut2simulate2
#$ -m abe
#$ -M ab08028
rundate=`date +%Y%m%d`
model=elut2.msmc.trimAncient.8.644114e-09
mkdir -p ${model}
for j in {1..6}
do
mkdir -p ${model}/group_$j.${model}
cd ${model}/group_$j.${model}
cp ../../../macs ./
cp ../../../msformatter ./
for i in {1..10}
do
# northern sea otter msmc model -- trimmed msmc
mu=8.644114e-09
r=1e-08
Na=12244.1746867
rho=0.000489766987469
theta=0.000423360167312
date=`date +%Y%m%d`
SEED=$((date+$RANDOM+((j-1)*10)+i))
# this is a new addition! need to have a different random seed for each simulation; if they start within a second of each other, they will have the same seed. not an issue for big simulations of 30Mb because those are slow, but 100kb can start within a second of each other!
./macs 2 30000000 -t 0.000423360167312 -r 0.000489766987469 -s $SEED -eN 0.0 0.31369842092 -eN 0.0076989765492 0.230281508202 -eN 0.0155979482953 0.198226326898 -eN 0.023707473624 0.214922772457 -eN 0.032039386431 0.27401235463 -eN 0.040606087505 0.35569903322 -eN 0.049420804354 0.436984654093 -eN 0.0584988903355 0.498223984644 -eN 0.0678564074235 0.530676045769 -eN 0.0775110710305 0.536368673347 -eN 0.0874822500075 0.510272725996 -eN 0.097791439055 0.510272725996 -eN 0.108462494928 0.455502738352 -eN 0.119521636438 0.455502738352 -eN 0.130998153067 0.404181175726 -eN 0.142924877378 0.404181175726 -eN 0.155338657431 0.367509179736 -eN 0.168280829187 0.367509179736 -eN 0.181798161335 0.346825490052 -eN 0.195944744936 0.346825490052 -eN 0.210781284802 0.342140865472 -eN 0.2263793512 0.342140865472 -eN 0.242821616055 0.355583907267 -eN 0.260201616745 0.355583907267 -eN 0.27863745602 0.393780841558 -eN 0.298263771015 0.393780841558 -eN 0.319243543525 0.473306970094 -eN 0.34177991028 0.473306970094 -eN 0.366120887055 0.636426341833 -eN 0.392580627165 0.636426341833 -eN 0.421563042015 1.0  > group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt 
./msformatter < group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt
done
cd ../../
done
