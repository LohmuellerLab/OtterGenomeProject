################## Making a new plotting script for sliding windows ##########################
require(ggplot2)
require(RColorBrewer)
cols <- RColorBrewer::brewer.pal(name="Set1",n=8)
colsrep <- rep(cols,220/8)
elutCol="#377EB8"
pbraCol="#4DAF4A"
elut2Col="#C77CFF"
window="1000000"
windowsize=as.numeric(as.character(window))
step="20000"
filter=0.8 # may need to change this

multiScaffHet_elut <- function(inputFilePath) {
  input <- read.table(inputFilePath,header=T,stringsAsFactors = F)
  input$hets_divBy_calls <- input$hets_GIDGET/input$calls_GIDGET
  return(input)
}

multiScaffHet_pbra <- function(inputFilePath) {
  input <- read.table(inputFilePath,header=T,stringsAsFactors = F)
  input$hets_divBy_calls <- input$hets_PteBra_1/input$calls_PteBra_1
  return(input)
}

# divides het by column in output that calls_PREFIX 
multiScaffHet_elut2 <- function(inputFilePath) {
  input <- read.table(inputFilePath,header=T,stringsAsFactors = F)
  input$hets_divBy_calls <- input$hets_Elfin/input$calls_Elfin
  return(input)
}

########### use .fai file to do chunk nums #####
fai <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/ferretGenome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa.fai")
head(fai)
colnames(fai) <- c("scaffold","length")
# get chunk numbers (row number in .fai file)
fai1_220 <- fai[1:220,]
fai1_220$chunk <- row.names(fai1_220) # note that after chunk 220 this is inaccurate 

#### Get position shift per chromosome: can use this as lookup table ######
posshift <- head(c(0,cumsum(fai1_220$length)),-1) # adds up chromsomes cumulatively ; cuts off the last one because you're going to shift stuff down
### that is what leads to the NAs. Why? I guess because that doesn't add it in?
# aha this step shifts everything down one so that the "0" box is the shift for hte first chromosome, etc.
names(posshift) <- as.factor(fai1_220$scaffold) # this has the first scaffold have a shift of 0, the next scaffold has a shift of length of the first scaffold, the third scaffold ahs a shift of len scaff 1+2, etc.

############### Read in data (50/250 DP Filter) ###########
# these are from totally filtered HQ sites vcf files
filename_elut="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/SlidingWindowheterozygosity/results_20180212_50_250DPFilter_1Mb/01_Elut_CA_Gidget.raw_variants.1-220Scaffs.20171006.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz_het_1000000win_20000step.txt"
filename_pbra="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/slidingWindowHeterozygosity/results_20170212_50_250DPFilter_1Mb/PteBra_1.raw_variants.1-220Scaffs.20171206.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz_het_1000000win_20000step.txt"
# this was final filtering using only lib 283
filename_elut2="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/genome-stats/slidingWindowHeterozygosity/02_Elut_SEAK_Elfin.raw_variants.1-220Scaffs.20190215_nso_lib283only.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz_het_1000000win_20000step.txt"

scaffs1_220_elut <- multiScaffHet_elut(filename_elut) 
scaffs1_220_pbra <- multiScaffHet_pbra(filename_pbra)
scaffs1_220_elut2 <- multiScaffHet_elut2(filename_elut2) 
######### Merge with .fai chunk info #########
# add chunk info:
scaffs1_220_chunks_elut <- merge(scaffs1_220_elut,fai1_220[,c("scaffold","chunk")],by.x="chromo",by.y="scaffold")
dim(scaffs1_220_chunks_elut) # 93000 windows
scaffs1_220_chunks_pbra <- merge(scaffs1_220_pbra,fai1_220[,c("scaffold","chunk")],by.x="chromo",by.y="scaffold")
dim(scaffs1_220_chunks_pbra) # 93006 windows
scaffs1_220_chunks_elut2 <- merge(scaffs1_220_elut2,fai1_220[,c("scaffold","chunk")],by.x="chromo",by.y="scaffold")
dim(scaffs1_220_chunks_elut2) # 93003 windows


############### RBIND elut and pbra #################
# make colnames the same:
colnames(scaffs1_220_chunks_pbra)[8:9] <- c("calls","hets")
colnames(scaffs1_220_chunks_elut)[8:9] <- c("calls","hets")
colnames(scaffs1_220_chunks_elut2)[8:9] <- c("calls","hets")
scaffs1_220_chunks_pbra$spp <- "Giant Otter"
scaffs1_220_chunks_elut$spp <- "S. Sea Otter"
scaffs1_220_chunks_elut2$spp <- "N. Sea Otter"
scaffs1_220_chunks <- rbind(scaffs1_220_chunks_elut,scaffs1_220_chunks_elut2,scaffs1_220_chunks_pbra) # both otters
dim(scaffs1_220_chunks)


################## Add filter flag to windows (instead of new dataframe) ###################
scaffs1_220_chunks$passFlag <- "fail" # don't use NA, creates weird behavior
scaffs1_220_chunks[scaffs1_220_chunks$calls >= (windowsize*filter),]$passFlag <- "pass"
# also get relative position for each window that is halfway between start and stop and is shifted by pos shift
scaffs1_220_chunks$genpos <- (scaffs1_220_chunks$window_start +0.5*(scaffs1_220_chunks$window_end - scaffs1_220_chunks$window_start)) + posshift[scaffs1_220_chunks$chromo]

############# add in mean heterozygosity as lines ###########
meanElutHet <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/overallHeterozygosity/elut.GenomeWideHet.50_250DPFilter.txt",sep=":")
meanPbraHet <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/overall-heterozygosity/pbra.GenomeWideHet.50_250DPFilter.txt",sep=":")
meanElut2Het <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/genome-stats/overall-heterozygosity/elut2.GenomeWideHet.50_250DPFilter.283.txt",sep=":")
############ Input data for plotting ###############
### windows passing filters
inputdata_pbra=scaffs1_220_chunks[scaffs1_220_chunks$spp=="Giant Otter" & scaffs1_220_chunks$chunk %in% seq(1,220) & scaffs1_220_chunks$passFlag=="pass",]
inputdata_elut=scaffs1_220_chunks[scaffs1_220_chunks$spp=="S. Sea Otter" & scaffs1_220_chunks$chunk %in% seq(1,220) & scaffs1_220_chunks$passFlag=="pass",]
inputdata_elut2=scaffs1_220_chunks[scaffs1_220_chunks$spp=="N. Sea Otter" & scaffs1_220_chunks$chunk %in% seq(1,220) & scaffs1_220_chunks$passFlag=="pass",]


faildata_pbra=scaffs1_220_chunks[scaffs1_220_chunks$spp=="Giant Otter" & scaffs1_220_chunks$chunk %in% seq(1,220) & scaffs1_220_chunks$passFlag=="fail",]
faildata_elut=scaffs1_220_chunks[scaffs1_220_chunks$spp=="S. Sea Otter" & scaffs1_220_chunks$chunk %in% seq(1,220) & scaffs1_220_chunks$passFlag=="fail",]
faildata_elut2=scaffs1_220_chunks[scaffs1_220_chunks$spp=="N. Sea Otter" & scaffs1_220_chunks$chunk %in% seq(1,220) & scaffs1_220_chunks$passFlag=="fail",]
############# Plot overlaid on each other: ######################

plot_relpos <- ggplot(data=inputdata_pbra,aes(x=genpos,y=hets_divBy_calls,color=spp))+
  geom_bar(size=0.1,alpha=1,stat="identity")+
  geom_bar(data=inputdata_elut,aes(x=genpos,y=hets_divBy_calls,color=spp),size=0.1,alpha=0.5,stat="identity")+
  geom_bar(data=inputdata_elut2,aes(x=genpos,y=hets_divBy_calls,color=spp),size=0.1,alpha=0.5,stat="identity")+
  xlab("Position Along Ferret Reference Genome (Gb)")+
  ylab("Heterozygotes/bp")+
  ggtitle("Sliding Window Heterozygosity\nScaffolds 1-220")+
  theme_classic() +
  scale_x_continuous(expand = c(0, 0),breaks=c(0.25e09,0.5e09,0.75e09,1e09,1.25e09,1.5e09,1.75e09,2e09),labels=c("0.25","0.5","0.75","1","1.25","1.5","1.75","2"))+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none")+
  scale_color_manual(values = c(pbraCol,elutCol,elut2Col))+
  scale_fill_manual(values = c(pbraCol,elut2Col,elutCol))+
  theme(legend.title=element_blank())+
  geom_hline(aes(yintercept=meanElutHet$V2),color="black")+
  geom_hline(aes(yintercept=meanPbraHet$V2),color="black",linetype="dashed")
### Be careful with labels if you change scale!! 
plot_relpos


############################# save plots ###########################
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/1Mb.Elut.PBra.220scaffs",as.character(filter),".0perc.long.lines.relativePosition.50_250DPFilter.pdf",sep=""),plot_relpos,height=5,width=9,device="pdf")
# also save as .png:
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/1Mb.Elut.PBra.220scaffs",as.character(filter),".0perc.long.lines.relativePosition.50_250DPFilter.png",sep=""),plot_relpos,height=3,width=9,device="png")
# also save as tif:
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/1Mb.Elut.PBra.220scaffs",as.character(filter),".0perc.long.lines.relativePosition.50_250DPFilter.tif",sep=""),plot_relpos,height=3,width=9,device="tiff",compression="xlm")
######################plotting each population separately in a stack: #############################
# need to add averages to the plot
## attempt to rbind?
allSpeciesSW <- rbind(inputdata_pbra,inputdata_elut,inputdata_elut2)
levels(as.factor(allSpeciesSW$spp))
# reorder levels so it goes pbra > elut1 > elut2
allSpeciesSW$spp <-  factor(allSpeciesSW$spp, levels = c("Giant Otter", "S. Sea Otter", "N. Sea Otter"))


# make auxiliary data frame of averages to add lines to each facet separately
means <- data.frame(spp=c("Giant Otter","S. Sea Otter","N. Sea Otter"),meanHet=c(meanPbraHet$V2,meanElutHet$V2,meanElut2Het$V2))
levels(allSpeciesSW$spp)
levels(allSpeciesSW$spp) <- as.factor(c("Giant Otter" , "S. Sea Otter" "N. Sea Otter"))
plotall3 <- ggplot(allSpeciesSW,aes(x=genpos,y=hets_divBy_calls,color=spp))+
  geom_bar(size=0.1,alpha=0.9,stat="identity")+
  #geom_hline(data=means, aes(yintercept = meanHet),linetype="dashed")+
  facet_wrap(~spp,ncol = 1)+
  xlab("Position Along Ferret Reference Genome (Gb)")+
  ylab("Heterozygotes/bp")+
  ggtitle("Sliding Window Heterozygosity\nScaffolds 1-220")+
  theme_bw() +
  scale_x_continuous(expand = c(0, 0),breaks=c(0.25e09,0.5e09,0.75e09,1e09,1.25e09,1.5e09,1.75e09,2e09),labels=c("0.25","0.5","0.75","1","1.25","1.5","1.75","2"))+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none")+
  scale_color_manual(values = c(pbraCol,elutCol,elut2Col))+
  scale_fill_manual(values = c(pbraCol,elutCol,elut2Col))+
  theme(legend.title=element_blank())
plotall3
# looks better without average

ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/1Mb.Elut1.Elut2.PBra.220scaffs",as.character(filter),".0perc.long.lines.relativePosition.correctDPFilters.noAvgLines.pdf",sep=""),plotall3,height=5,width=10,device="pdf",dpi=300)
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/1Mb.Elut1.Elut2.PBra.220scaffs",as.character(filter),".0perc.long.lines.relativePosition.correctDPFilters.noAvgLines.tif",sep=""),plotall3,height=5,width=10,device="tiff",dpi=300)
# also save as .png:
ggsave(paste("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/1Mb.Elut.PBra.220scaffs",as.character(filter),".0perc.long.lines.relativePosition.correctDPFilters.noAvgLines.png",sep=""),plotall3,height=5,width=10,device="png",dpi=300)



