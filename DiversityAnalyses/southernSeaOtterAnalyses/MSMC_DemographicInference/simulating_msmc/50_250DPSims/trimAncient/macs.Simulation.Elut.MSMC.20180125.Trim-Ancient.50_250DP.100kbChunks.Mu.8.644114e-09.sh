#!/bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=28G,highp
#$ -N elutMSMCTrim
#$ -m abe
#$ -M ab08028
rundate=`date +%Y%m%d`
model=elut.msmc.trimAncient.100kbChunks.8.644114e-09
mkdir -p ${model}
for j in {1..18}
do
mkdir -p ${model}/group_$j.${model}
cd ${model}/group_$j.${model}
cp ../../../macs ./
cp ../../../msformatter ./
for i in {1..1000}
do
# sea otter msmc model -- full msmc
mu=8.644114e-09
r=1e-08
Na=22759.3264133
rho=0.000910373056531
theta=0.000786936848318
date=`date +%Y%m%d`
SEED=$((date+$RANDOM+((j-1)*1000)+i))
# this is a new addition! need to have a different random seed for each simulation; if they start within a second of each other, they will have the same seed. not an issue for big simulations of 30Mb because those are slow, but 100kb can start within a second of each other!
./macs 2 100000 -t 0.000786936848318 -r 0.000910373056531 -s $SEED -eN 0.0 0.088103054401 -eN 0.003617749005 0.0339730033926 -eN 0.0073294699725 0.0929355795355 -eN 0.011140182365 0.14648499415 -eN 0.0150553377 0.185969764821 -eN 0.019080692475 0.218631178707 -eN 0.023222829175 0.244857652103 -eN 0.02748860985 0.262357053343 -eN 0.031885659 0.270286941254 -eN 0.036422363575 0.269710116013 -eN 0.041107872975 0.257097104104 -eN 0.04595209905 0.257097104104 -eN 0.05096647855 0.232676303912 -eN 0.0561630836 0.232676303912 -eN 0.061556019525 0.208900140554 -eN 0.067160281175 0.208900140554 -eN 0.072993531975 0.191032772099 -eN 0.07907496025 0.191032772099 -eN 0.085426804125 0.18040432147 -eN 0.09207422445 0.18040432147 -eN 0.09904606725 0.177336636081 -eN 0.1063754991 0.177336636081 -eN 0.11410140495 0.182994441404 -eN 0.122268896425 0.182994441404 -eN 0.13093172625 0.200959926622 -eN 0.140153559 0.200959926622 -eN 0.1500120375 0.240224202955 -eN 0.160602468 0.240224202955 -eN 0.172039218 0.323417652254 -eN 0.18447350675 0.323417652254 -eN 0.1980921345 0.512531459732 -eN 0.2131479805 0.512531459732 -eN 0.2299777935 1.0  > group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt 
./msformatter < group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt
done
cd ../../
done
