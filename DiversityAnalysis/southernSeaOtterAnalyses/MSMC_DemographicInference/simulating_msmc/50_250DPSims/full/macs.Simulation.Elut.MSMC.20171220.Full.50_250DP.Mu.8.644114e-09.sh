#!/bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=28G,highp
#$ -N elutMSMC1
#$ -m abe
#$ -M ab08028
rundate=`date +%Y%m%d`
model=elut.msmc.full.8.644114e-09
mkdir -p ${model}
for j in {1..6}
do
mkdir -p ${model}/group_$j.${model}
cd ${model}/group_$j.${model}
cp ../../../macs ./
cp ../../../msformatter ./
for i in {1..10}
do
# sea otter msmc model -- full msmc
mu=8.644114e-09
r=1e-08
Na=411054.933124
rho=0.016442197325
theta=0.0142128228087
date=`date +%Y%m%d`
SEED=$((date+$RANDOM+((j-1)*10)+i))
# this is a new addition! need to have a different random seed for each simulation; if they start within a second of each other, they will have the same seed. not an issue for big simulations of 30Mb because those are slow, but 100kb can start within a second of each other!
./macs 2 30000000 -t 0.0142128228087 -r 0.016442197325 -s $SEED -eN 0.0 0.00487809781987 -eN 0.00020030785146 0.00188102029959 -eN 0.00040581875097 0.00514566550505 -eN 0.00061681061658 0.00811059429737 -eN 0.0008335852884 0.0102967906221 -eN 0.0010564614927 0.0121051907162 -eN 0.0012858036891 0.013557300448 -eN 0.0015219918162 0.0145262088657 -eN 0.001765448028 0.014965271611 -eN 0.0020166366939 0.0149333338993 -eN 0.0022760643987 0.0142349755244 -eN 0.0025442799426 0.0142349755244 -eN 0.0028219165566 0.0128828424686 -eN 0.0031096426512 0.0128828424686 -eN 0.0034082392113 0.0115664017228 -eN 0.0037185364731 0.0115664017228 -eN 0.0040415124267 0.0105771196633 -eN 0.004378229493 0.0105771196633 -eN 0.0047299189545 0.00998864265535 -eN 0.0050979739194 0.00998864265535 -eN 0.005483991537 0.00981879077556 -eN 0.0058898081772 0.00981879077556 -eN 0.0063175768254 0.0101320526483 -eN 0.0067697952261 0.0101320526483 -eN 0.007249439565 0.0111267672455 -eN 0.007760034828 0.0111267672455 -eN 0.00830587995 0.0133007552199 -eN 0.008892251856 0.0133007552199 -eN 0.009525482856 0.017907017584 -eN 0.010213945671 0.017907017584 -eN 0.010967983074 0.028377887842 -eN 0.011801596506 0.028377887842 -eN 0.012733431102 0.0553680897108 -eN 0.013789871487 0.0553680897108 -eN 0.015009474393 0.13764305417 -eN 0.016451974611 0.13764305417 -eN 0.018217422639 0.44745963203 -eN 0.02049346593 0.44745963203 -eN 0.023701414176 1.0 -eN 0.029185405713 1.0  > group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt 
./msformatter < group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt
done
cd ../../
done
