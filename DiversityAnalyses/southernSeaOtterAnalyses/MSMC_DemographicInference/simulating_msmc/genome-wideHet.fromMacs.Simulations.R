#### GENOME WIDE HETEROZYGOSITY FROM SIMULATIONS
# this input is the number of segregating sites
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/DemographicInference_WholeGenome/msmc")
require(ggplot2)
windowsize=30000000 # 30Mb windows in this simulation

# model 1 is full msmc trajectory
# model 2 is with ancient removed
# model 3 is with modern and ancient removed
input1 <- read.table("simulatingMSMC/simulations_20171031_1.25mu/elut.msmc.full.20171031.SegSiteCounts.txt",header=F)
input2 <- read.table("simulatingMSMC/simulations_20171031_1.25mu/elut.msmc.trimAncient.20171031.SegSiteCounts.txt",header=F)
input3 <- read.table("simulatingMSMC/simulations_20171031_1.25mu/elut.msmc.trimModernAndAncient.20171031.SegSiteCounts.txt",header=F)
windowCount <- dim(input1)[1]
### genome-wide het 1:
het1= sum(input1)/(windowCount*windowsize) # 0.0008306267
### genome-wide het 2:
het2= sum(input2)/(windowCount*windowsize) # 9.895905e-05
### genome-wide het 3:
het3 = sum(input3)/(windowCount*windowsize) # 9.751095e-05

## from the data:
hetdata=0.0002

hets <- rbind(het1,het2,het3,hetdata)
hets <- data.frame(hets)
hets$model <- row.names(hets)
hets$label <- NA
hets[hets$model=="het1",]$label <- "full trajectory"
hets[hets$model=="het2",]$label <- "trimmed anc."
hets[hets$model=="het3",]$label <- "trim anc. + modern"
hets[hets$model=="hetdata",]$label <- "data"

## 
ggplot(hets,aes(y=hets,x=label))+
  geom_bar(stat="identity")
