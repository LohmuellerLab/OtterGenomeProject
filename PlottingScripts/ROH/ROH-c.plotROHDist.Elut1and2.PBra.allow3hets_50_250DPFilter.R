# Want to plot results:
# catted together in step roh-b
#setwd("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/runsOfHomozygosity")
# 20190130 adding NSO
require(ggplot2)
elutCol="#377EB8"
pbraCol="#4DAF4A"
elut2Col="#C77CFF"
## ** use the 220 ** 
### UPDATED 20180216 to use 50/250 DP Results
#rohResults_elut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/runsOfHomozygosity/results_20171025/scaffs1_220.catted.hom",header=T,strip.white =T)
rohResults_elut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/runsOfHomozygosity/results_20180216_50_250DPFilter/01_Elut_CA_Gidget.scaffs1_220.catted.hom",header=T,strip.white =T)
#rohResults_pbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/runsOfHomozygosity/scaffs1_220.catted.hom",header=T,strip.white =T)
rohResults_pbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/runsOfHomozygosity/results_20180216_50_250DPFilter/PteBra_1.scaffs1_220.catted.hom",header=T,strip.white =T)
rohResults_elut2 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/genome-stats/runsOfHomozygosity/02_Elut_SEAK_Elfin.scaffs1_220.catted.hom",header=T,strip.white =T) # this uses new 20x-250x dp filter 20180201
dim(rohResults_pbra)
dim(rohResults_elut)
dim(rohResults_elut2)
# adding spaces for better plotting
rohResults_pbra$spp <- " Giant Otter"
rohResults_elut$spp <- " S. Sea Otter"
rohResults_elut2$spp <- " N. Sea Otter"
rohResults <- rbind(rohResults_pbra,rohResults_elut,rohResults_elut2)
rohResults$MB <- as.numeric(rohResults$KB) / 1000 # to go from KB to MB 
head(rohResults)
p <- ggplot(rohResults,aes(x=MB,color=spp,fill=spp))+
  geom_histogram(binwidth=1,breaks=seq(0,max(round(rohResults$MB)),by=2),size=0.5,position="dodge")+
  #breaks=seq(0,15,by=1)
  stat_bin(breaks=seq(0,max(round(rohResults$MB)),by=2),aes(y=..count.., label=..count..), geom="text", vjust=-.2,position="dodge") +
  xlab("ROH Size (Mb)")+
  ggtitle("Runs of Homozygosity (ROH)\nScaffolds 1-220\n3 Hets. allowed per window")+
  theme_bw()+
  scale_color_manual(values=c(pbraCol,elut2Col,elutCol))+
  scale_fill_manual(values=c(pbraCol,elut2Col,elutCol))+
  scale_y_continuous()+ # so that label shows up 
  scale_x_continuous(breaks=seq(0,max(round(rohResults$MB)),by=2),labels=seq(0,max(round(rohResults$MB)),by=2))+
  theme(axis.text = element_text(size=14))
p

# want to show total amount of MB in each category. Hm, how to do that?
# make 3 categories:
rohResults$MBCategory <- NA
rohResults[rohResults$MB < 5,]$MBCategory <- "< 5Mb"
rohResults[rohResults$MB > 10,]$MBCategory <- "> 10Mb"
rohResults[rohResults$MB <= 10 & rohResults$MB >= 5,]$MBCategory <- "5-10Mb"
library(dplyr)
library(tidyr)
rohMBTally <- rohResults %>% 
  group_by(MBCategory,spp) %>% 
  summarise(MB = sum(MB)) %>%
  ungroup %>%
  complete(spp,MBCategory,fill = list(MB = 0))
rohMBTally
rohMBTally <- data.frame(rohMBTally)
# reorder levels for plotting
levels(as.factor(rohMBTally$MBCategory))
rohMBTally$MBCategory <-  factor(rohMBTally$MBCategory, levels = c("< 5Mb" ,"5-10Mb", "> 10Mb"))
levels(rohMBTally$MBCategory)

levels(as.factor(rohMBTally$spp))
rohMBTally$spp <-  factor(rohMBTally$spp, levels = c(" Giant Otter" ," S. Sea Otter", " N. Sea Otter"))
levels(rohMBTally$spp)

write.table(rohMBTally,"/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/ROH/rohMBTally.50_250DPFilter.includesNSO.txt",row.names=F,quote=F)

p2 <- ggplot(rohMBTally,aes(x=MBCategory,y=MB,color=spp,fill=spp))+
  theme_bw()+
  geom_bar(stat="identity",position="dodge")+
  xlab("ROH Size (Mb)")+
  theme(legend.title = element_blank())+
  ylab("Mb contained in ROHs")+
  ggtitle("Runs of Homozygosity (ROH)")+
  scale_color_manual(values=c(pbraCol,elutCol,elut2Col))+
  scale_fill_manual(values=c(pbraCol,elutCol,elut2Col))+
  theme(axis.text = element_text(size=14),legend.text = element_text(size=14))+
  theme(legend.position= c(0.7,0.8),legend.background = element_rect(fill="transparent"))
  
p2

ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/ROH/roh.plot.1-220.3hetsAllowed.ElutPbra.50_250DPFilter.newElut2Filter.includesNSO.pdf",p2,width=4,height=5,device="pdf",units="in")

