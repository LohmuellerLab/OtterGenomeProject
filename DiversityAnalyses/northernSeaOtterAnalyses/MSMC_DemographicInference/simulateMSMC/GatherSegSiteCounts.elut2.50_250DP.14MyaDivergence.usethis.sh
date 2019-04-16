# can run in the shell
### this has been modified to account for zero site simulations

spp=elut2
wd=/u/flashscratch/a/ab08028/msmcSimulations/50_250DPSimulations_14MyaDivergence/$spp

reasonableMu=8.644114e-09
#20180828 update: reasonable mu based on 14.4 Mya divergence time from 5-cal phylogeny (4 year gen time)

mu=$reasonableMu
mods="full trimAncient"


mkdir -p $wd/segSiteCountsAllModels/

for mod in $mods
do
model=${spp}.msmc.${mod}.${mu}
echo $model
echo "group block segsites" > $wd/segSiteCountsAllModels/${model}.SegSiteCounts.txt
for j in {1..6}
do
for i in {1..10}
do
# if you can grep seg sites, do so
file=$wd/$mod/$model/group_${j}.${model}/group_${j}_block_${i}.${model}.msFormat.OutputFile.*.txt
if grep -q segsites $file
then
grep segsites $file | awk -v j=$j -v i=$i '{OFS=" ";print j,i,$2}' >> $wd/segSiteCountsAllModels/${model}.SegSiteCounts.txt
else
echo "$j $i 0" >> $wd/segSiteCountsAllModels/${model}.SegSiteCounts.txt
fi
done
done
done

