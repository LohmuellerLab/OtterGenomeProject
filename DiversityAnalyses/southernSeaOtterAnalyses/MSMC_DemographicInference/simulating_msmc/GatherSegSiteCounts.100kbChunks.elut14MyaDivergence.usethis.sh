#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=2G
#$ -N sims
#$ -m abe
#$ -M ab08028
### this has been modified to account for zero site simulations
# important to note! if no segsites, it was zero sites.
reasonableMu=8.644114e-09
#20180828 update: reasonable mu based on 14.4 Mya divergence time from 5-cal phylogeny (4 year gen time)

spp="elut"
mod1="msmc.trimAncient.rmBneck.100kbChunks"
mod2="msmc.trimAncient.100kbChunks"
wd=/u/flashscratch/a/ab08028/msmcSimulations/50_250DPSimulations_14MyaDivergence/$spp


mkdir -p $wd/segSiteCountsAllModels/

# trim and remove bneck:
mu=$reasonableMu
model1=${spp}.${mod1}.${mu}
echo $model1
echo "group block segsites" > $wd/segSiteCountsAllModels/$model1.SegSiteCounts.txt
for j in {1..18}
do
for i in {1..1000}
do
echo $j ":" $i 
# if you can grep seg sites, do so
file=$wd/trimAncient.rmBneck/$model1/group_${j}.${model1}/group_${j}_block_${i}.${model1}.msFormat.OutputFile.*.txt
if grep -q segsites $file
then
grep segsites $file | awk -v j=$j -v i=$i '{OFS=" ";print j,i,$2}' >> $wd/segSiteCountsAllModels/$model1.SegSiteCounts.txt
else
echo "$j $i 0" >> $wd/segSiteCountsAllModels/$model1.SegSiteCounts.txt
fi
done
done


# trim and keep bneck:
model2=${spp}.${mod2}.${mu}
echo $model2
echo "group block segsites" > $wd/segSiteCountsAllModels/$model2.SegSiteCounts.txt
for j in {1..18}
do
for i in {1..1000}
do
echo $j ":" $i 
# if you can grep seg sites, do so
file=$wd/trimAncient/$model2/group_${j}.${model2}/group_${j}_block_${i}.${model2}.msFormat.OutputFile.*.txt
if grep -q segsites $file
then
grep segsites $file | awk -v j=$j -v i=$i '{OFS=" ";print j,i,$2}' >> $wd/segSiteCountsAllModels/$model2.SegSiteCounts.txt
else
echo "$j $i 0" >> $wd/segSiteCountsAllModels/$model2.SegSiteCounts.txt
fi
done
done
