# 20180828: adding statistical test
######## Step 7_C 
require(reshape2)
require(ggplot2)
require(scales)
require(ggpubr)
elutCol="#377EB8"
pbraCol="#4DAF4A"
elut2Col="#C77CFF"
fullpagewidth_mm=170
halfpagewidth_mm=85
maxheight=225
################ EMPIRICAL VALUES ##############
########### vep results ############
########## Calculating: het/hom of deleterious variants #######
elut_inputCounts <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/includesShared/01_Elut_CA_Gidget.20171006.VepResultsCountSummary.txt",header=T) 
misHetElut <- elut_inputCounts[elut_inputCounts$category=="missense" & elut_inputCounts$genotype=="Heterozygous",]$count # missense het count
# make a dataframe of sites with that count and category (derived allele count is 1 for hets)
misHomElut <- elut_inputCounts[elut_inputCounts$category=="missense" & elut_inputCounts$genotype=="HomozygousAlt",]$count
synHetElut <- elut_inputCounts[elut_inputCounts$category=="synonymous" & elut_inputCounts$genotype=="Heterozygous",]$count
synHomElut <- elut_inputCounts[elut_inputCounts$category=="synonymous" & elut_inputCounts$genotype=="HomozygousAlt",]$count
sgHetElut <- elut_inputCounts[elut_inputCounts$category=="stopgained" & elut_inputCounts$genotype=="Heterozygous",]$count
sgHomElut <- elut_inputCounts[elut_inputCounts$category=="stopgained" & elut_inputCounts$genotype=="HomozygousAlt",]$count


pbra_inputCounts <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/includesShared/loadCalcs_includesSharedVars/PteBra_1.20171206.VepResultsCountSummary.IncludesShared.txt",header=T) 
misHetPbra <- pbra_inputCounts[pbra_inputCounts$category=="missense" & pbra_inputCounts$genotype=="Heterozygous",]$count # missense het count
# make a dataframe of sites with that count and category (derived allele count is 1 for hets)
misHomPbra <- pbra_inputCounts[pbra_inputCounts$category=="missense" & pbra_inputCounts$genotype=="HomozygousAlt",]$count
synHetPbra <- pbra_inputCounts[pbra_inputCounts$category=="synonymous" & pbra_inputCounts$genotype=="Heterozygous",]$count
synHomPbra <- pbra_inputCounts[pbra_inputCounts$category=="synonymous" & pbra_inputCounts$genotype=="HomozygousAlt",]$count
sgHetPbra <- pbra_inputCounts[pbra_inputCounts$category=="stopgained" & pbra_inputCounts$genotype=="Heterozygous",]$count
sgHomPbra <- pbra_inputCounts[pbra_inputCounts$category=="stopgained" & pbra_inputCounts$genotype=="HomozygousAlt",]$count


# NSO:

elut2_inputCounts <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/VEP_50_250DPFilter/VEP-50_250DPFilter_canon_cds/includesShared/loadCalcs_includesSharedVars/filteredVepOutput2-ForLoadCalcs/02_Elut_SEAK_Elfin.20190215_nso_lib283only.VepResultsCountSummary.txt",header=T) 
misHetElut2 <- elut2_inputCounts[elut2_inputCounts$category=="missense" & elut2_inputCounts$genotype=="Heterozygous",]$count # missense het count
# make a dataframe of sites with that count and category (derived allele count is 1 for hets)
misHomElut2 <- elut2_inputCounts[elut2_inputCounts$category=="missense" & elut2_inputCounts$genotype=="HomozygousAlt",]$count
synHetElut2 <- elut2_inputCounts[elut2_inputCounts$category=="synonymous" & elut2_inputCounts$genotype=="Heterozygous",]$count
synHomElut2 <- elut2_inputCounts[elut2_inputCounts$category=="synonymous" & elut2_inputCounts$genotype=="HomozygousAlt",]$count
sgHetElut2 <- elut2_inputCounts[elut2_inputCounts$category=="stopgained" & elut2_inputCounts$genotype=="Heterozygous",]$count
sgHomElut2 <- elut2_inputCounts[elut2_inputCounts$category=="stopgained" & elut2_inputCounts$genotype=="HomozygousAlt",]$count

# Totals ?: Hom + 1/2 Het
totalMisPbra = (as.numeric(misHomPbra[1])+0.5*as.numeric(misHetPbra[1]))
totalSynPbra =  (as.numeric(synHomPbra[1])+0.5*as.numeric(synHetPbra[1]))
totalMisElut = (as.numeric(misHomElut[1])+0.5*as.numeric(misHetElut[1]))
totalSynElut =  (as.numeric(synHomElut[1])+0.5*as.numeric(synHetElut[1]))
totalMisElut2 = (as.numeric(misHomElut2[1])+0.5*as.numeric(misHetElut2[1]))
totalSynElut2 =  (as.numeric(synHomElut2[1])+0.5*as.numeric(synHetElut2[1]))
############# total coding sites #########
elut_codingCallableSites_perChunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/calculatingCodingCallableSites/countCallableCDSRegions_50_250DPFilter//elut.callableSites.ferretCDSRegions.txt",header=T)
total_Coding_Callable_elut <- sum(elut_codingCallableSites_perChunk$callableCDSSites) 
total_Coding_Callable_elut # 30680418 (formerly 25,591,577 (callable coding sites) (~10% of genome, good.) )
pbra_codingCallableSites_perChunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/calculatingCodingCallableSites/countCallableCDSRegions_50_250DPFilter/pbra.callableSites.ferretCDSRegions.txt",header=T)
total_Coding_Callable_pbra <- sum(pbra_codingCallableSites_perChunk$callableCDSSites) 
total_Coding_Callable_pbra # 30070013 (forermly 29,200,510 AHA! Pbra had more coding callable sites.)

elut2_codingCallableSites_perChunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret//genome-stats/CallableCDSSites/elut2.callableSites.ferretCDSRegions.txt",header=T)
total_Coding_Callable_elut2 <- sum(elut2_codingCallableSites_perChunk$callableCDSSites) 
total_Coding_Callable_elut2 # 31002191 * more sites because higher coverage * *this is with improved 20x-250x Dp filter *


average_Coding_Callable <- (total_Coding_Callable_elut+total_Coding_Callable_pbra+total_Coding_Callable_elut2)/3
average_Coding_Callable

############ Synonymous:  ############ 
# add syn. het sites plus 2* syn hom sites together
# divide by num callable sites to get a proportion 
e1_prop <- (as.numeric(synHetElut[1])+2*as.numeric(synHomElut[1]))/(total_Coding_Callable_elut)
p1_prop <-  (as.numeric(synHetPbra[1])+2*as.numeric(synHomPbra[1]))/(total_Coding_Callable_pbra)
# nso: 
e1_prop2 <- (as.numeric(synHetElut2[1])+2*as.numeric(synHomElut2[1]))/(total_Coding_Callable_elut2)
e1_prop # gidget
p1_prop
e1_prop2 # nso
# then rescale 
SynDerivedElutRescaled <- e1_prop*average_Coding_Callable
SynDerivedPbraRescaled <- p1_prop*average_Coding_Callable
SynDerivedElut2Rescaled <- e1_prop2*average_Coding_Callable

SynDerivedElutRescaled
SynDerivedPbraRescaled
SynDerivedElut2Rescaled
############ Missense (Non-Synonymous) :  ############ 
e2_prop <- (as.numeric(misHetElut[1])+2*as.numeric(misHomElut[1]))/total_Coding_Callable_elut
p2_prop <-  (as.numeric(misHetPbra[1])+2*as.numeric(misHomPbra[1]))/total_Coding_Callable_pbra
# nso:
e2_prop2 <- (as.numeric(misHetElut2[1])+2*as.numeric(misHomElut2[1]))/total_Coding_Callable_elut2

e2_prop
p2_prop
e2_prop2

# then rescale 
MisDerivedElutRescaled <- e2_prop*average_Coding_Callable
MisDerivedPbraRescaled <- p2_prop*average_Coding_Callable
MisDerivedElut2Rescaled <- e2_prop2*average_Coding_Callable

MisDerivedElutRescaled
MisDerivedPbraRescaled
MisDerivedElut2Rescaled
############ LOF #######
e3_prop <- (as.numeric(sgHetElut[1])+2*as.numeric(sgHomElut[1]))/total_Coding_Callable_elut
p3_prop <-  (as.numeric(sgHetPbra[1])+2*as.numeric(sgHomPbra[1]))/total_Coding_Callable_pbra
e3_prop2 <- (as.numeric(sgHetElut2[1])+2*as.numeric(sgHomElut2[1]))/total_Coding_Callable_elut2

e3_prop
p3_prop
e3_prop2
SGDerivedElutRescaled <- e3_prop*average_Coding_Callable
SGDerivedPbraRescaled <- p3_prop*average_Coding_Callable
SGDerivedElut2Rescaled <- e3_prop2*average_Coding_Callable

SGDerivedElutRescaled
SGDerivedPbraRescaled
SGDerivedElut2Rescaled
########### Make a dataframe: #######
dfDerived <- data.frame(rbind(SGDerivedElutRescaled,SGDerivedElut2Rescaled,SGDerivedPbraRescaled,MisDerivedElutRescaled,MisDerivedElut2Rescaled,MisDerivedPbraRescaled,SynDerivedElutRescaled,SynDerivedElut2Rescaled,SynDerivedPbraRescaled))
colnames(dfDerived) <- "RescaledCount"
dfDerived$spp <- c("S. Sea Otter","N. Sea Otter","Giant Otter","S. Sea Otter","N. Sea Otter","Giant Otter","S. Sea Otter","N. Sea Otter","Giant Otter")
dfDerived$type <- c("Stop-Gained","Stop-Gained","Stop-Gained","Missense","Missense","Missense","Synonymous","Synonymous","Synonymous") 
write.table(dfDerived,"/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/geneticLoad/CountsOfRescaledDerivedAlleles.byType.txt",sep="\t",quote=F,row.names=F)
# plot just the estimates before boots
ggplot(dfDerived,aes(x=spp,y=RescaledCount,fill=spp))+
  geom_bar(stat="identity")+
  facet_wrap(~type,scales="free")
###################### BOOTSTRAPS ########################
# these are new bootstraps as of 20190201 that include the average based on all three otters
elutBoots <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/includesShared/elut_derivedAllele_bootstraps/elut.allBootstrapsConcatted.results.txt",header=T) # got from Hoffman, result of step 7b
dim(elutBoots) # 1000 bootstraps
# melt it:

elutBoots_melt <- melt(elutBoots)
elutBoots_melt$spp <- "S. Sea Otter"
elutBoots_melt$type <- NA
elutBoots_melt[elutBoots_melt$variable=="MissenseDerived",]$type <- "Missense"
elutBoots_melt[elutBoots_melt$variable=="SynonymousDerived",]$type <- "Synonymous"
elutBoots_melt[elutBoots_melt$variable=="Stop_GainedDerived",]$type <- "Stop-Gained"
# derived allele counts only:
colnames(elutBoots_melt) <- c("variable","RescaledCount","spp","type")
# combine with original:
elutBoots_melt_derived <- elutBoots_melt[elutBoots_melt$variable %in% c("MissenseDerived","SynonymousDerived","Stop_GainedDerived"),]

####### pbra:
pbraBoots <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/includesShared/loadCalcs_includesSharedVars/pbra_derivedAllele_bootstraps/pbra.allBootstrapsConcatted.results.txt",header=T) # got from Hoffman, result of step 7b THESE INCLUDE SHARED ALLELEs
dim(pbraBoots) # 1000 bootstraps
# melt it:

pbraBoots_melt <- melt(pbraBoots)
pbraBoots_melt$spp <- "Giant Otter"
pbraBoots_melt$type <- NA
pbraBoots_melt[pbraBoots_melt$variable=="MissenseDerived",]$type <- "Missense"
pbraBoots_melt[pbraBoots_melt$variable=="SynonymousDerived",]$type <- "Synonymous"
pbraBoots_melt[pbraBoots_melt$variable=="Stop_GainedDerived",]$type <- "Stop-Gained"
# derived allele counts only:
colnames(pbraBoots_melt) <- c("variable","RescaledCount","spp","type")
# combine with original:
pbraBoots_melt_derived <- pbraBoots_melt[pbraBoots_melt$variable %in% c("MissenseDerived","SynonymousDerived","Stop_GainedDerived"),]


####### nothern sea otter  ##########################
elut2Boots <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/VEP_50_250DPFilter/VEP-50_250DPFilter_canon_cds/includesShared/loadCalcs_includesSharedVars/filteredVepOutput2-ForLoadCalcs/elut2_derivedAllele_bootstraps/elut2.allBootstrapsConcatted.results.txt",header=T) # got from Hoffman, result of step 7b THESE INCLUDE SHARED ALLELEs
dim(elut2Boots) # 1000 bootstraps
# melt it:

elut2Boots_melt <- melt(elut2Boots)
elut2Boots_melt$spp <- "N. Sea Otter"
elut2Boots_melt$type <- NA
elut2Boots_melt[elut2Boots_melt$variable=="MissenseDerived",]$type <- "Missense"
elut2Boots_melt[elut2Boots_melt$variable=="SynonymousDerived",]$type <- "Synonymous"
elut2Boots_melt[elut2Boots_melt$variable=="Stop_GainedDerived",]$type <- "Stop-Gained"
# derived allele counts only:
colnames(elut2Boots_melt) <- c("variable","RescaledCount","spp","type")
# combine with original:
elut2Boots_melt_derived <- elut2Boots_melt[elut2Boots_melt$variable %in% c("MissenseDerived","SynonymousDerived","Stop_GainedDerived"),]

########### rbind sea otter and giant otter bootstraps ######
allsppBoots_melt_derived <- rbind(elutBoots_melt_derived,elut2Boots_melt_derived,pbraBoots_melt_derived)

######## rearrange factors: #######
print(levels(as.factor(allsppBoots_melt_derived$type)))
allsppBoots_melt_derived$type <- factor(allsppBoots_melt_derived$type,levels=c("Synonymous","Missense","Stop-Gained"))
print(levels(as.factor(allsppBoots_melt_derived$type)))

print(levels(as.factor(dfDerived$type)))
dfDerived$type <- factor(dfDerived$type,levels=c("Synonymous","Missense","Stop-Gained"))
print(levels(as.factor(dfDerived$type)))

############## I. Plot derived allele count by category (for manuscript) #############
# reorder levels to match other plots (doesn't reassign, just changes order of plotting)
levels(as.factor(allsppBoots_melt_derived$spp))
allsppBoots_melt_derived$spp <- factor(allsppBoots_melt_derived$spp,levels=c("Giant Otter" , "S. Sea Otter" ,"N. Sea Otter"))
levels(as.factor(allsppBoots_melt_derived$spp))
p1 <- ggplot(allsppBoots_melt_derived,aes(x=spp,y=RescaledCount/1000,fill=spp))+
  theme_bw()+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values=c(pbraCol,elutCol,elut2Col))+
  geom_point(data=dfDerived,aes(x=spp,y=RescaledCount/1000),color="black",size=2)+
  facet_wrap(~type,scales = "free_y")+
  theme(strip.text.x = element_text(size = 10, colour = "black"),legend.title = element_blank(),legend.position = "none",axis.text.x=element_text(size=6))+
  ylab("Derived Alleles\n(x1000)")+
  xlab("")
p1


############## II. Plot the number of homozygous derived GENOTYPES (rescaled) for manuscript ###########
# Rescale homozygous genotype counts:
########## synonymous:
SynHomDerivedGTsElutRescaled <- (synHomElut/total_Coding_Callable_elut)*average_Coding_Callable
########## missense:
MisHomDerivedGTsElutRescaled <- (misHomElut/total_Coding_Callable_elut)*average_Coding_Callable
########## stop-gained
SGHomDerivedGTsElutRescaled <- (sgHomElut/total_Coding_Callable_elut)*average_Coding_Callable

########## synonymous:
SynHomDerivedGTsPbraRescaled <- (synHomPbra/total_Coding_Callable_pbra)*average_Coding_Callable
########## missense:
MisHomDerivedGTsPbraRescaled <- (misHomPbra/total_Coding_Callable_pbra)*average_Coding_Callable
########## stop-gained
SGHomDerivedGTsPbraRescaled <- (sgHomPbra/total_Coding_Callable_pbra)*average_Coding_Callable


########## synonymous:
SynHomDerivedGTsElut2Rescaled <- (synHomElut2/total_Coding_Callable_elut2)*average_Coding_Callable
########## missense:
MisHomDerivedGTsElut2Rescaled <- (misHomElut2/total_Coding_Callable_elut2)*average_Coding_Callable
########## stop-gained
SGHomDerivedGTsElut2Rescaled <- (sgHomElut2/total_Coding_Callable_elut2)*average_Coding_Callable

######### make a hom derived dataframe #######

dfHomDerivedGTs <- data.frame(rbind(SGHomDerivedGTsElutRescaled,SGHomDerivedGTsElut2Rescaled,SGHomDerivedGTsPbraRescaled,MisHomDerivedGTsElutRescaled,MisHomDerivedGTsElut2Rescaled,MisHomDerivedGTsPbraRescaled,SynHomDerivedGTsElutRescaled,SynHomDerivedGTsElut2Rescaled,SynHomDerivedGTsPbraRescaled))

colnames(dfHomDerivedGTs) <- "RescaledCount"
dfHomDerivedGTs$spp <- c("S. Sea Otter","N. Sea Otter","Giant Otter","S. Sea Otter","N. Sea Otter","Giant Otter","S. Sea Otter","N. Sea Otter","Giant Otter")
dfHomDerivedGTs$type <- c("Stop-Gained","Stop-Gained","Stop-Gained","Missense","Missense","Missense","Synonymous","Synonymous","Synonymous") 
write.table(dfHomDerivedGTs,"/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/geneticLoad/HomozygousDerived.Alleles.ByCategory.txt", row.names=F,quote=F,sep="\t")
############## Hom Derived GTs bootstraps ##########
# 20180220 kirk said to just plot homozygous derived  genotypes across 3 categories
pbraBoots_melt[pbraBoots_melt$variable=="missenseHomCount",]$type <- "Missense"
pbraBoots_melt[pbraBoots_melt$variable=="synHomGenotypeCount",]$type <- "Synonymous"
pbraBoots_melt[pbraBoots_melt$variable=="sgHomCount",]$type <- "Stop-Gained"
# just homozygous genotypes:
pbraBoots_melt_HomGTs <- pbraBoots_melt[pbraBoots_melt$variable %in% c("missenseHomCount","synHomGenotypeCount","sgHomCount"),]

elutBoots_melt[elutBoots_melt$variable=="missenseHomCount",]$type <- "Missense"
elutBoots_melt[elutBoots_melt$variable=="synHomGenotypeCount",]$type <- "Synonymous"
elutBoots_melt[elutBoots_melt$variable=="sgHomCount",]$type <- "Stop-Gained"
elutBoots_melt_HomGTs <- elutBoots_melt[elutBoots_melt$variable %in% c("missenseHomCount","synHomGenotypeCount","sgHomCount"),]
## note that the bootstrap counts are already rescaled because the avg number of sites was drawn with replacement

# 
elut2Boots_melt[elut2Boots_melt$variable=="missenseHomCount",]$type <- "Missense"
elut2Boots_melt[elut2Boots_melt$variable=="synHomGenotypeCount",]$type <- "Synonymous"
elut2Boots_melt[elut2Boots_melt$variable=="sgHomCount",]$type <- "Stop-Gained"
elut2Boots_melt_HomGTs <- elut2Boots_melt[elut2Boots_melt$variable %in% c("missenseHomCount","synHomGenotypeCount","sgHomCount"),]

# all spp just keeping it named both spp for convenience 
allsppBoots_melt_HomGTs <- rbind(elutBoots_melt_HomGTs,pbraBoots_melt_HomGTs,elut2Boots_melt_HomGTs)

##############  rearrange factors: 
print(levels(as.factor(allsppBoots_melt_HomGTs$type)))
allsppBoots_melt_HomGTs$type <- factor(allsppBoots_melt_HomGTs$type,levels=c("Synonymous","Missense","Stop-Gained"))
print(levels(as.factor(allsppBoots_melt_HomGTs$type)))

print(levels(as.factor(dfHomDerivedGTs$type)))
dfHomDerivedGTs$type <- factor(dfHomDerivedGTs$type,levels=c("Synonymous","Missense","Stop-Gained"))
print(levels(as.factor(dfHomDerivedGTs$type)))


########### plot ###########
levels(as.factor(allsppBoots_melt_HomGTs$spp))
allsppBoots_melt_HomGTs$spp <- factor(allsppBoots_melt_HomGTs$spp,levels=c("Giant Otter" , "S. Sea Otter" ,"N. Sea Otter"))
levels(as.factor(allsppBoots_melt_HomGTs$spp))

p2 <- ggplot(allsppBoots_melt_HomGTs,aes(x=spp,y=RescaledCount/1000,fill=spp))+
  theme_bw()+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values=c(pbraCol,elutCol,elut2Col))+
  geom_point(data=dfHomDerivedGTs,aes(x=spp,y=RescaledCount/1000),color="black",size=2)+
  facet_wrap(~type,scales = "free_y")+
  theme(strip.text.x = element_text(size = 10, colour = "black"),legend.title = element_blank(),legend.position = "none",axis.text.x=element_text(size=6))+
  ylab("Homozygous Derived\nGenotypes (x1000)")+
  xlab("")
p2
# expect elut higher across all categories
# maybe changed stopgained to stopgained in a domain?

############### arrange plots in a grid ###########
require(gridExtra)
grid.arrange(p1,p2)
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/VEP-GeneticLoad/loadFigure.derived.homozygous.INCLUDESNSO.pdf",grid.arrange(p1,p2),device="pdf",width=16.9,height=8,units = "cm",dpi=300)
# save a tiff:
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/VEP-GeneticLoad/loadFigure.derived.homozygous.INCLUDESNSO.tif",grid.arrange(p1,p2),device="tiff",width=16.9,height=8,units = "cm",dpi=300,compression = "lzw")
# save as a tiff file 
############################## Z TEST FOR SIGNIFICANCE #############
######### Significance testing ##########
# Calculate a z score for bootstrap
# z = [(point estimate (elut) - point estimate (pbra) )- (truemean1 - truemean2)]/sqrt(variance (from bootstrap) 1+variance (from bootstrap) 2) 
# (truemean1 - truemean2) <-- this part gets set to 0 (assume true mean1=truemean2)
# and you don't divide variance by n. 
require(dplyr)
#################### *** YOU ARE HERE *** ######################
# want to use point estimates instead of mean; don't get mean from bootstraps
summary_allsppBoots_melt_derived <- allsppBoots_melt_derived %>%
  group_by(spp,type) %>%
  summarise(variance=var(RescaledCount))
# okay then want to get the z Score for each one
#df=data.frame(category=,zScore=integer(),pval2sided=integer())
mylist=list()
for(i in c("Synonymous","Missense","Stop-Gained")){
   #print(i)
   varDF1 <- summary_allsppBoots_melt_derived %>%
    filter(spp=="Giant Otter" & type==i)
   varDF2 <- summary_allsppBoots_melt_derived %>%
     filter(spp=="S. Sea Otter" & type==i)
   varDF3 <- summary_allsppBoots_melt_derived %>%
     filter(spp=="N. Sea Otter" & type==i)
   # get point estimate from data (rescaled)
   ptEst1 <- dfDerived %>%
     filter(spp=="Giant Otter" & type==i)
   ptEst2 <- dfDerived %>%
     filter(spp=="S. Sea Otter" & type==i)
   ptEst3 <- dfDerived %>%
     filter(spp=="N. Sea Otter" & type==i)
   # s. sea otter point estimate - giant otter point estimate
   xdiff1=ptEst2$RescaledCount - ptEst1$RescaledCount
   # s. sea otter point estimate - n. sea otter point estimate
   xdiff2=ptEst2$RescaledCount - ptEst3$RescaledCount
   # n. sea otter point estimate - giant otter point estimate
   xdiff3=ptEst3$RescaledCount - ptEst1$RescaledCount
   sqvarSum1=sqrt(varDF1$variance + varDF2$variance)
   sqvarSum2=sqrt(varDF2$variance + varDF3$variance)
   sqvarSum3=sqrt(varDF3$variance + varDF1$variance)
   zScore1=xdiff1/sqvarSum1
   zScore2=xdiff2/sqvarSum2
   zScore3=xdiff3/sqvarSum3
   print(paste(i," Z-Score1 (giant otter & s. sea otter): ",zScore1,sep=""))
   print(paste(i," Z-Score2 (s. sea otter & n. sea otter): ",zScore2,sep=""))
   print(paste(i," Z-Score3 (giant otter & n. sea otter): ",zScore3,sep=""))
   pval2sided1=2*pnorm(-abs(zScore1))
   pval2sided2=2*pnorm(-abs(zScore2))
   pval2sided3=2*pnorm(-abs(zScore3))
# put into a vector
   vec <- c(i,zScore1,pval2sided1,zScore2,pval2sided2,zScore3,pval2sided3)
   print(paste(i," two-tailed p-value1 (giant otter & s. sea otter): ",pval2sided1,sep=""))
   print(paste(i," two-tailed p-value2 (s. sea otter & n. sea otter): ",pval2sided2,sep=""))
   print(paste(i," two-tailed p-value3 (giant otter & n. sea otter): ",pval2sided3,sep=""))
   mylist[[i]] <- vec
}
df <- do.call("rbind",mylist)
colnames(df) <- c("category","Z-Score1 (giant otter & s. sea otter)","Two-Sided_Pvalue1 (giant otter & s. sea otter)","Z-Score2 (s. sea otter & n. sea otter)","Two-Sided_Pvalue2 (s. sea otter & n. sea otter)","Z-Score3 (giant otter & n. sea otter)","Two-Sided_Pvalue3 (giant otter & n. sea otter)")
write.table(df,"/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/geneticLoad/ZTest.GeneticLoad.DerivedAlleles.txt",row.names=F,quote=F,sep="\t")
# note you have to ensure Z score is negative. If it's negative then it is giving you the probability of being less than that score because lower.tail=T is the default (means prob of being less than) so you use -abs() to ensure it's negative then you
# multiply by 2 to get two-tailed because I am just hypothesizing a difference, not a direction
### Do for genotypes:

summary_allsppBoots_melt_HomGTs <- allsppBoots_melt_HomGTs %>%
  group_by(spp,type) %>%
  summarise(variance=var(RescaledCount))

mylist2=list()

for(i in c("Synonymous","Missense","Stop-Gained")){
  #print(i)
  varDF1 <- summary_allsppBoots_melt_HomGTs %>%
    filter(spp=="Giant Otter" & type==i)
  varDF2 <- summary_allsppBoots_melt_HomGTs %>%
    filter(spp=="S. Sea Otter" & type==i)
  varDF3 <- summary_allsppBoots_melt_HomGTs %>%
    filter(spp=="N. Sea Otter" & type==i)
  # get point estimate from data (rescaled)
  ptEst1 <- dfHomDerivedGTs %>%
    filter(spp=="Giant Otter" & type==i)
  ptEst2 <- dfHomDerivedGTs %>%
    filter(spp=="S. Sea Otter" & type==i)
  ptEst3 <- dfHomDerivedGTs %>%
    filter(spp=="N. Sea Otter" & type==i)
  # s. sea otter point estimate - giant otter point estimate
  xdiff1=ptEst2$RescaledCount - ptEst1$RescaledCount
  # s. sea otter point estimate - n. sea otter point estimate
  xdiff2=ptEst2$RescaledCount - ptEst3$RescaledCount
  # n. sea otter point estimate - giant otter point estimate
  xdiff3=ptEst3$RescaledCount - ptEst1$RescaledCount
  sqvarSum1=sqrt(varDF1$variance + varDF2$variance)
  sqvarSum2=sqrt(varDF2$variance + varDF3$variance)
  sqvarSum3=sqrt(varDF3$variance + varDF1$variance)
  zScore1=xdiff1/sqvarSum1
  zScore2=xdiff2/sqvarSum2
  zScore3=xdiff3/sqvarSum3
  print(paste(i," Z-Score1 (giant otter & s. sea otter): ",zScore1,sep=""))
  print(paste(i," Z-Score2 (s. sea otter & n. sea otter): ",zScore2,sep=""))
  print(paste(i," Z-Score3 (giant otter & n. sea otter): ",zScore3,sep=""))
  pval2sided1=2*pnorm(-abs(zScore1))
  pval2sided2=2*pnorm(-abs(zScore2))
  pval2sided3=2*pnorm(-abs(zScore3))
  # put into a vector
  vec <- c(i,zScore1,pval2sided1,zScore2,pval2sided2,zScore3,pval2sided3)
  print(paste(i," two-tailed p-value1 (giant otter & s. sea otter): ",pval2sided1,sep=""))
  print(paste(i," two-tailed p-value2 (s. sea otter & n. sea otter): ",pval2sided2,sep=""))
  print(paste(i," two-tailed p-value3 (giant otter & n. sea otter): ",pval2sided3,sep=""))
  mylist2[[i]] <- vec
}
df2 <- do.call("rbind",mylist2)
colnames(df2) <-  c("category","Z-Score1 (giant otter & s. sea otter)","Two-Sided_Pvalue1 (giant otter & s. sea otter)","Z-Score2 (s. sea otter & n. sea otter)","Two-Sided_Pvalue2 (s. sea otter & n. sea otter)","Z-Score3 (giant otter & n. sea otter)","Two-Sided_Pvalue3 (giant otter & n. sea otter)")



write.table(df2,"/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/geneticLoad/ZTest.GeneticLoad.DerivedGenotypes.txt",row.names=F,quote=F,sep="\t")

