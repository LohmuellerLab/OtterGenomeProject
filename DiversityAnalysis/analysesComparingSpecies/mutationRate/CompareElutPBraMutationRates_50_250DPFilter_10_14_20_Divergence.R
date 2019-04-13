# Compare mutation rates of elut and pbra:
grid_elut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/mutationRate/elut.GridOfMutationRates.50_250DPFilter.10_14_20_MyaDivergence.txt",header=T)
grid_pbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/mutationRate/pbra.GridOfMutationRates.50_250DPFilter.10_14_20_MyaDivergence.txt",header=T)
grid_elut$spp <- "Sea Otter"
grid_pbra$spp <- "Giant Otter"

grids <- rbind(grid_elut,grid_pbra)


require(ggplot2)
ggplot(grids,aes(x=divergence_yrs,y=mutationRate,color=as.factor(generation_time),shape=spp))+
  theme_bw()+
  geom_point(size=3)
  
# My choices: average between elut and pbra (they are very close)
high <- mean(grids[grids$divergence_yrs== 10400000 & grids$generation_time==10,]$mutationRate)
high
# High: 9mya + 10yr/gen (extreme high mutation rate) = average(3.558263e-08, 3.542660e-08)
# Reasonable:  11mya, 4yr/gen = average()
##### 20180813: now I have the phylogeny. Divergence of otters from ferret estimate is 12.7 ; low = 9.3, high = 16.9. So I'm making the 3 values 9, 13, 17 (keep low and high the same; changing reasonable)
#reasonable <- mean(grids[grids$divergence_yrs== 11e06 & grids$generation_time==4,]$mutationRate)
reasonable <- mean(grids[grids$divergence_yrs== 14400000 & grids$generation_time==4,]$mutationRate)
reasonable
# Low : 17 mya, 2yr/gen
low <- mean(grids[grids$divergence_yrs== 20100000 & grids$generation_time==2,]$mutationRate)
low

# write these out. and make combo grid somehow
high 
reasonable
low

write.table(paste("high_mutation_rate ",high,"\nreasonable_mutation_rate ",reasonable,"\nlow_mutation_rate ",low,sep=""),"/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/MutationRateGrid/averageMutationRateGrid.high.reasonable.low.10_14_20_MyaDivergence.txt",row.names = F,col.names=F,quote=F)
