#### Calculate mutation rate 
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/")
# final variants:
variants <- read.table("overallHeterozygosity/step-by-step-counts_50_250DPFilter/counts.8.Vars.5.finalPASS")
# final hets:
hets <- read.table("overallHeterozygosity/step-by-step-counts_50_250DPFilter/counts.8.Vars.5.finalPASS.hets")
homs <- read.table("overallHeterozygosity/step-by-step-counts_50_250DPFilter/counts.8.Vars.5.finalPASS.homs")
# final sites (snps, biallelic, indels removed, bad sites removed etc)
totalsitesinput <- read.table("overallHeterozygosity/step-by-step-counts_50_250DPFilter/counts.9.HQSites.rmBadVars")
# this includes homozygous ancestral, heterozygous, homozygous derived
homDerived <- sum(homs$V1)
mutations <- homDerived + sum(hets)/2
mutations
mutations_bp <- mutations/sum(totalsitesinput)
### Divergence times in years
#### 20180829: updating this based on the 5-calibration phylogeny (independent clock)
# 95% credible intervals:

time1=10.4e6 # lower bound (10.4)

time2=14.4e6 

time3=20.1e6

times=c(time1,time2,time3)
### possible generation times
gen1 = 4 # years/gen
gen2 = 2
gen3= 8
gen4=10
gen5=6 # cite gagne 2018
# Ralls uses 4, 8 10 as sea otter gneration times. Also putting in two in case ancestor had much slower generations


gens=c(gen1,gen2,gen3,gen4,gen5)
### make a grid:
grid <- expand.grid(x=times,y=gens)
colnames(grid) <- c("divergence_yrs","generation_time")
### get divergence in generations:
grid$divergence_gen <- grid$divergence_yrs / grid$generation_time
grid$mutationRate <- mutations_bp / grid$divergence_gen

write.table(grid,"mutationRate/elut.GridOfMutationRates.50_250DPFilter.10_14_20_MyaDivergence.txt",row.names=F,quote=F)


# calculate:
min(grid$mutationRate)
max(grid$mutationRate)
mean(grid$mutationRate)
require(ggplot2)
ggplot(grid,aes(y=mutationRate,x=divergence_yrs,color=as.factor(generation_time)))+
  geom_point()
