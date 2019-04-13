setwd("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp")
require(ggplot2)
require(ggrepel)
require(RColorBrewer)
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
require(qvalue)
require(clipr)
########### Colors #######
#display.brewer.pal(name="Set1",n=6)
palette <- brewer.pal(name="Set1",n=6)
artifactcol = palette[1]
realcol = palette[3]
unexaminedcol <- "slategrey"
filtercol1 <- palette[4]
filtercol2 <- palette[2]
################## GENE INFO ##################
### Need: Elut gene info (can be used for any set of alignments):
allElutGenes_fullInfo <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Annotation/MakerRound5_blast1e-06_USETHIS_Nov11/FULLINFO_GradedGeneSets_Nov17/mRNA_FullInfo_withGrades.Nov17.txt",sep=";",header=T,quote="",comment.char="")
# add a jbrowse column:
allElutGenes_fullInfo$jbrowse <- paste(allElutGenes_fullInfo$scaffold,":",allElutGenes_fullInfo$start,"-",allElutGenes_fullInfo$stop,sep="")
library(stringr)
head(allElutGenes_fullInfo)
int1e <- str_split_fixed(allElutGenes_fullInfo$blastInfo,":",2)
head(int1e)
allElutGenes_fullInfo$similarToSymbol <- str_split_fixed(int1e[,1],"to ",2)[,2]
allElutGenes_fullInfo$similarToSymbol

### in step 8b (only once) I got gene symbols for all the genes in the cfam and hsap lookup tables
# so this can be used for any orthologs, because all had to come from those sets
# I got the symbols from gene convert and manually added them using Google Sheets (do not use Excel because of name conversion issues!)

cfamtt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T)
hsaptt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T,quote="")

############################## FUNCTIONS ######################################
# adding to this function to get 3 reps and deal with convergence
codemlGetPVals_allInfo_guidance_3Reps <- function(date,flag,humanSymbols,dogSymbols){
  require(qvalue)
  # read in file: (marked by flag and date), omit NA:
  codemldf.1 <- na.omit(read.table(paste("output_codeml_",date,"/codeml_",flag,".rep.1.",date,".fullInfo.txt",sep=""),header=T))
  codemldf.2 <- na.omit(read.table(paste("output_codeml_",date,"/codeml_",flag,".rep.2.",date,".fullInfo.txt",sep=""),header=T))
  codemldf.3 <- na.omit(read.table(paste("output_codeml_",date,"/codeml_",flag,".rep.3.",date,".fullInfo.txt",sep=""),header=T))
  # omit rows with NA (paml didn't finish right); but keep those with symbol as N/A (just missing symbol)
  # Deal with 3 replicates before doing the filters below:
  # merge the three reps in 2 steps:
  TwoReps <- merge(codemldf.1,codemldf.2,by="TName",suffixes = c("",".rep2"))
  ThreeReps <- merge(TwoReps,codemldf.3,by="TName",suffixes = c("",".rep3"))
  # get max of each one:
  headerNull=(paste(flag,".lnL_null",sep=""))
  headerAlt=(paste(flag,".lnL_alt",sep=""))
  max_null <- apply(ThreeReps[,c(headerNull,paste(headerNull,"rep2",sep="."),paste(headerNull,"rep3",sep="."))],1,FUN=max)
  max_alt <- apply(ThreeReps[,c(headerAlt,paste(headerAlt,"rep2",sep="."),paste(headerAlt,"rep3",sep="."))],1,FUN=max)
  ThreeReps_MAXlnls <- cbind(ThreeReps,max_null,max_alt) # this gets the max values attached
  # calculate maxLRT:
  ThreeReps_MAXlnls$maxLRT <- 2*(ThreeReps_MAXlnls$max_alt - ThreeReps_MAXlnls$max_null)
  # Set remainder to zero:
  ThreeReps_MAXlnls[ThreeReps_MAXlnls$maxLRT < 0,]$maxLRT <- 0
  
  ########NOW DO EVERYTHING ELSE USING maxLRT:
  # 20170925: adding 120bp length filter
  ThreeReps_MAXlnls <- ThreeReps_MAXlnls[ThreeReps_MAXlnls$seqLength_mayInclGaps > 120,]
  # 20170925: optional dramatic move, but get rid of alignments with >1 significant BEB sites (indicative of bad alignment)
  ThreeReps_MAXlnls <- ThreeReps_MAXlnls[ThreeReps_MAXlnls[,paste(flag,".bebCount",sep="")] <=1,]
  # convert to P value: df = 1
  ThreeReps_MAXlnls$uncorrPvalues <- (1-pchisq(ThreeReps_MAXlnls[,"maxLRT"],1))
  # want cluster info (names of other spp transcripts)
  # using CAT for this run because FNR was off in awk. In future use hsap.
  ThreeReps_MAXlnls_dogClusters <- merge(ThreeReps_MAXlnls,unique(dogSymbols[,c("cfam","Symbol")]),by.x="TName",by.y="cfam",all=F)
  ThreeReps_MAXlnls_humanClusters <- merge(ThreeReps_MAXlnls,unique(humanSymbols[,c("Transcript","Symbol")]),by.x="TName",by.y="Transcript",all=F)
  ThreeReps_MAXlnls_allClusters <- rbind(ThreeReps_MAXlnls_dogClusters,ThreeReps_MAXlnls_humanClusters)
  # this also gets all other spp and the gene symbol
  # want to get jbrowse coordinates
  return(ThreeReps_MAXlnls_allClusters)
}
# then put in the full info clusters with pvalues for each one:
codemlGetQVals_allInfo_guidance_3Reps <- function(otters, elut, pbra,fdr=0.05){
  # for each spp. want to combine pvalues to do qvalue correction across all gene/branches
  allP <- rbind(otters[,c("branch","TName","uncorrPvalues")],elut[,c("branch","TName","uncorrPvalues")],pbra[,c("branch","TName","uncorrPvalues")])
  qobj <- qvalue(allP$uncorrPvalues, pi0.meth="bootstrap", fdr.level=fdr)
  allP$qvalues  <- qobj$qvalues
  allP$pvaluesFromqvalue_check <- qobj$pvalues
  allP$significant <- qobj$significant
  return(allP)
}
# this is if you just want q value for one branch
# input is the result of codemlGetPVals_allInfo_guidance_3Reps()
codemlGetQVals_1branch <- function(input,fdr=0.05){
  qobj <- qvalue(input$uncorrPvalues, pi0.meth="bootstrap", fdr.level=fdr)
  input$qvalues  <- qobj$qvalues
  input$pvaluesFromqvalue_check <- qobj$pvalues
  input$significant <- qobj$significant
  return(input)
}
# this will give you info on each branch+TName
# then you can merge it with the others
gatherSignificantResults_scriptMaker <- function(sigResults,branch,date,note){
  rep="1.2.3" # all 3 reps
  inspect <- sigResults # whatever subset of the genes you want to inspect 
  # fix group names and cluster names by adding leading zeros
  inspect$group <- sprintf("%02d",inspect$group)
  inspect$cluster <- sprintf("%05s",inspect$cluster)
  # get rid of slashes 
  inspect$Symbol <- gsub("N/A","N_A",inspect$Symbol)
  # set up empty file so you don't just keep appending:
  write(paste("#",Sys.Date(),"This will gather the following files for visualization"),file=paste("output_codeml_",date,"/inspect_codeml_subset_",date,"_",branch,"_rep",rep,"_",note,".sh",sep=""))
  write(paste("mkdir -p $SCRATCH/visualInspection/",branch,"/","codeml_output_",date,"_",note,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"_",inspect$Symbol,sep=""),file=paste("output_codeml_",date,"/inspect_codeml_subset_",date,"_",branch,"_rep",rep,"_",note,".sh",sep=""),append=T)
  # cp beb info 
  write(paste("cp -r group_",inspect$group,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"/codeml_branchSite_results_rep_*_",date,"/bebInfo/*",branch,"* $SCRATCH/visualInspection/",branch,"/","codeml_output_",date,"_",note,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"_",inspect$Symbol,sep=""),file=paste("output_codeml_",date,"/inspect_codeml_subset_",date,"_",branch,"_rep",rep,"_",note,".sh",sep=""),append=T)
  # cp gblocksmed html
  write(paste("cp -r group_",inspect$group,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"/guidance*/*fastagbmed.htm $SCRATCH/visualInspection/",branch,"/","codeml_output_",date,"_",note,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"_",inspect$Symbol,sep=""),file=paste("output_codeml_",date,"/inspect_codeml_subset_",date,"_",branch,"_rep",rep,"_",note,".sh",sep=""),append=T)
}




############## Filtering Method A (#5): post stringent swamp (cleandata=1) 20170717 ###############
# (filtering method A: swamp + cleandata=1, elut pbra otters) 
otters_20170717 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170717,flag="otters",humanSymbols = hsaptt,dogSymbols = cfamtt)
elut_20170717 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170717,flag="elut",humanSymbols = hsaptt,dogSymbols = cfamtt)
pbra_20170717 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170717,flag="pbra",humanSymbols = hsaptt,dogSymbols = cfamtt)


# do q value correction across all branches/genes:
qVals_20170717 <- codemlGetQVals_allInfo_guidance_3Reps(otters_20170717,elut_20170717,pbra_20170717,0.05)

# can then look at each by branch:
otters_20170717_q <- merge(otters_20170717,qVals_20170717[qVals_20170717$branch=="otters",],by="TName")
elut_20170717_q <- merge(elut_20170717,qVals_20170717[qVals_20170717$branch=="elut",],by="TName")
pbra_20170717_q <- merge(pbra_20170717,qVals_20170717[qVals_20170717$branch=="pbra",],by="TName")

######### Filtering Method B (#6):20170921 (filtering method b: gblocks + cleandata=1, only otters ) ########

otters_20170921 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170921,flag="otters",humanSymbols = hsaptt,dogSymbols = cfamtt)
# do q value correction across all genes (only 1 branch)
otters_20170921_q <- codemlGetQVals_1branch(otters_20170921)

###### Filtering Method C (#7):20170922 (filtering method C: gblocks + cleandata=0, elut pbra otters ) ######
otters_20170922 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170922,flag="otters",humanSymbols = hsaptt,dogSymbols = cfamtt)
elut_20170922 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170922,flag="elut",humanSymbols = hsaptt,dogSymbols = cfamtt)
pbra_20170922 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170922,flag="pbra",humanSymbols = hsaptt,dogSymbols = cfamtt)

# do q value correction across all branches/genes:
qVals_20170922 <- codemlGetQVals_allInfo_guidance_3Reps(otters_20170922,elut_20170922,pbra_20170922,0.05)

# can then look at each by branch:
otters_20170922_q <- merge(otters_20170922,qVals_20170922[qVals_20170922$branch=="otters",],by="TName")
elut_20170922_q <- merge(elut_20170922,qVals_20170922[qVals_20170922$branch=="elut",],by="TName")
pbra_20170922_q <- merge(pbra_20170922,qVals_20170922[qVals_20170922$branch=="pbra",],by="TName")

####################### RESULTS OF VISUAL INSPECTION ###########################
# for each branch and filter, inspect the genes w q < 0.1, and with p < 0.01 q > 0.1
####### vi: otters 20170921 #####
artifacts_otters_20170921 <- read.table("output_codeml_20170921/inspection_otters1/artifacts_otters_20170921.txt",header=F)
colnames(artifacts_otters_20170921) <- "Symbol"

putativelyOkGenes_otters_20170921 <- read.table("output_codeml_20170921/inspection_otters1/putatively_good_otters_20170921.txt",header=F)
colnames(putativelyOkGenes_otters_20170921) <- "Symbol"

# label the genes:
otters_20170921_q$label <- "unexamined"
otters_20170921_q[otters_20170921_q$Symbol %in% putativelyOkGenes_otters_20170921$Symbol,]$label <- "putatively real"
otters_20170921_q[otters_20170921_q$Symbol %in% artifacts_otters_20170921$Symbol,]$label <- "artifact"


################### vi: otters 20170922 ###############
# otters 20170922
putativelyOkGenes_otters_20170922 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/putativelyReal_20170922.txt",header=F)
colnames(putativelyOkGenes_otters_20170922) <- "Symbol"
artifacts_otters_20170922 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/artifacts.20170922.txt",header=F)
colnames(artifacts_otters_20170922) <- "Symbol"
# label the genes:
otters_20170922_q$label <- "unexamined"
otters_20170922_q[otters_20170922_q$Symbol %in% putativelyOkGenes_otters_20170922$Symbol,]$label <- "putatively real"
otters_20170922_q[otters_20170922_q$Symbol %in% artifacts_otters_20170922$Symbol,]$label <- "artifact"

################### vi: otters 20170717 ###############
putativelyOkGenes_otters_20170717 <- read.table("output_codeml_20170717/otters_inspection1/putativelyReal_20170717.txt",header=F)
colnames(putativelyOkGenes_otters_20170717) <- "Symbol"
artifacts_otters_20170717 <- read.table("output_codeml_20170717/otters_inspection1/artifacts.20170717.txt",header=F)
colnames(artifacts_otters_20170717) <- "Symbol"
# label the genes:
otters_20170717_q$label <- "unexamined"
otters_20170717_q[otters_20170717_q$Symbol %in% putativelyOkGenes_otters_20170717$Symbol,]$label <- "putatively real"
otters_20170717_q[otters_20170717_q$Symbol %in% artifacts_otters_20170717$Symbol,]$label <- "artifact"

################### vi: elut 20170922 ###############
# elut
putativelyOKGenes_elut_20170922 <- read.table("output_codeml_20170922/inspection_elut/elut.20170922.3reps.putativelyGOOD.txt",header=T)
artifacts_elut_20170922 <- read.table("output_codeml_20170922/inspection_elut/elut.20170922.3reps.ARTIFACT.txt",header=T)
# label the genes:
elut_20170922_q$label <- "unexamined"
elut_20170922_q[elut_20170922_q$Symbol %in% putativelyOKGenes_elut_20170922$Symbol,]$label <- "putatively real"
elut_20170922_q[elut_20170922_q$Symbol %in% artifacts_elut_20170922$Symbol,]$label <- "artifact"


################### vi: elut 20170717 ###############
# need to inspect the elut genes that are with p < 0.01 (both sig and non sig)
# for now only inspected the genes that weren't already inspected in 20170922's inspection (keeping real/artifact calls from 20170922)
putativelyOKGenes_elut_20170717 <- read.table("output_codeml_20170717/inspection_elut/elut.20170717.3reps.putativelyGOOD.txt",header=T) # note that these are only genes that aren't already inspected by 20170922
artifacts_elut_20170717 <- read.table("output_codeml_20170717/inspection_elut/elut.20170717.3reps.ARTIFACT.txt",header=T)
# label the genes:
elut_20170717_q$label <- "unexamined"
elut_20170717_q[elut_20170717_q$Symbol %in% putativelyOKGenes_elut_20170717$Symbol,]$label <- "putatively real"
elut_20170717_q[elut_20170717_q$Symbol %in% artifacts_elut_20170717$Symbol,]$label <- "artifact"

################### vi: pbra 20170922 ###############
putativelyOkGenes_pbra_20170922 <- read.table("output_codeml_20170922/inspection_pbra/putatively_good_pbra_20170922.txt",header=F)
artifacts_pbra_20170922 <- read.table("output_codeml_20170922/inspection_pbra/artifacts_pbra_20170922.txt",header=F)
colnames(artifacts_pbra_20170922) <- "Symbol"
colnames(putativelyOkGenes_pbra_20170922) <- "Symbol"
# label the genes:
pbra_20170922_q$label <- "unexamined"
pbra_20170922_q[pbra_20170922_q$Symbol %in% putativelyOkGenes_pbra_20170922$Symbol,]$label <- "putatively real"
pbra_20170922_q[pbra_20170922_q$Symbol %in% artifacts_pbra_20170922$Symbol,]$label <- "artifact"

################### vi: pbra 20170717 ###############
putativelyOkGenes_pbra_20170717 <- read.table("output_codeml_20170717/inspection_pbra/putatively_good_pbra_20170717.txt",header=F)
artifacts_pbra_20170717 <- read.table("output_codeml_20170717/inspection_pbra/artifacts_pbra_20170717.txt",header=F)
colnames(artifacts_pbra_20170717) <- "Symbol"
colnames(putativelyOkGenes_pbra_20170717) <- "Symbol"
# label the genes:
pbra_20170717_q$label <- "unexamined"
pbra_20170717_q[pbra_20170717_q$Symbol %in% putativelyOkGenes_pbra_20170717$Symbol,]$label <- "putatively real"
pbra_20170717_q[pbra_20170717_q$Symbol %in% artifacts_pbra_20170717$Symbol,]$label <- "artifact"


########################### PLOTS FOR MANUSCRIPT #######################
##################### Single Filter PLOTS ########################################
#### 20170922 plot 1: elut (gets rid of genes with symbol NA or with LRT <1 (for ease of plotting)) 
elut_20170922_plot <- ggplot(elut_20170922_q[elut_20170922_q$maxLRT>=1 & elut_20170922_q$Symbol!="N/A",],aes(x=Symbol,y=maxLRT,color=label))+
  geom_point()+
  theme_bw()+
  theme(legend.title = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1,size=5))+
  geom_hline(yintercept = min(elut_20170922_q[elut_20170922_q$qvalues < 0.1,]$maxLRT),color="blue")+
  geom_hline(yintercept = min(elut_20170922_q[elut_20170922_q$pvaluesFromqvalue_check<0.01 & elut_20170922_q$qvalues > 0.1,]$maxLRT),color="blue",linetype=2) +
  ggtitle("Sea Otter\n Gblocks Filter")+
  geom_text_repel(
    data = elut_20170922_q[elut_20170922_q$Symbol %in% putativelyOKGenes_elut_20170922$Symbol,],
    aes(label = Symbol),
    size = 3,
    box.padding = 0.25,
    point.padding = 0.3,
    color="black"
  )

elut_20170922_plot

#### 20170922 plot 2: pbra
pbra_20170922_plot <- ggplot(pbra_20170922_q[pbra_20170922_q$maxLRT>=1 & pbra_20170922_q$Symbol!="N/A",],aes(x=Symbol,y=maxLRT))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5))+
  geom_hline(yintercept = min(pbra_20170922_q[pbra_20170922_q$qvalues < 0.1,]$maxLRT),color="blue")+
  geom_hline(yintercept = min(pbra_20170922_q[pbra_20170922_q$pvaluesFromqvalue_check<0.01 & pbra_20170922_q$qvalues > 0.1,]$maxLRT),color="blue",linetype=2) 

pbra_20170922_plot

#### 20170922 plot 3: otters
otters_20170922_plot <- ggplot(otters_20170922_q[otters_20170922_q$maxLRT>=1 & otters_20170922_q$Symbol!="N/A",],aes(x=Symbol,y=maxLRT))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5))+
  geom_hline(yintercept = min(otters_20170922_q[otters_20170922_q$qvalues < 0.1,]$maxLRT),color="blue")+
  geom_hline(yintercept = min(otters_20170922_q[otters_20170922_q$pvaluesFromqvalue_check<0.01 & otters_20170922_q$qvalues > 0.1,]$maxLRT),color="blue",linetype=2) +
  geom_point(data=otters_20170922_q[otters_20170922_q$Symbol %in% putativelyOkGenes_otters7$V1,],aes(x=Symbol,y=(maxLRT)),color="green")+
  geom_point(data=otters_20170922_q[otters_20170922_q$Symbol %in% artifacts_otters7$V1,],aes(x=Symbol,y=(maxLRT)),color="red")

otters_20170922_plot

###################### COMPARE filter PLOTS ######################
####### A vs C (3 spp) #####
# OTTERS
AvsC_otters <- merge(otters_20170717_q,otters_20170922_q,by="TName",suffixes=c(".A",".C"))
# merge labels:
AvsC_otters$mergedLabel <- AvsC_otters$label.C
AvsC_otters[AvsC_otters$label.A=="artifact",]$mergedLabel <- "artifact"
AvsC_otters[AvsC_otters$label.A=="putatively real",]$mergedLabel <- "putatively real"

# ELUT
AvsC_elut <- merge(elut_20170717_q,elut_20170922_q,by="TName",suffixes=c(".A",".C"))
# merge labels:
AvsC_elut$mergedLabel <- AvsC_elut$label.C
AvsC_elut[AvsC_elut$label.A=="artifact",]$mergedLabel <- "artifact"
AvsC_elut[AvsC_elut$label.A=="putatively real",]$mergedLabel <- "putatively real"

# PBRA
AvsC_pbra <- merge(pbra_20170717_q,pbra_20170922_q,by="TName",suffixes=c(".A",".C"))
# merge labels:
AvsC_pbra$mergedLabel <- AvsC_pbra$label.C
AvsC_pbra[AvsC_pbra$label.A=="artifact",]$mergedLabel <- "artifact"
AvsC_pbra[AvsC_pbra$label.A=="putatively real",]$mergedLabel <- "putatively real"

# this is the minimum LRT value of genes in the category q < 0.1, or p < 0.01
minqA_elut <- min(AvsC_elut[AvsC_elut$qvalues.A <= 0.1,]$maxLRT.A)
minpA_elut <- min(AvsC_elut[AvsC_elut$qvalues.A > 0.1 & AvsC_elut$pvaluesFromqvalue_check.A < 0.01,]$maxLRT.A)
minqC_elut <- min(AvsC_elut[AvsC_elut$qvalues.C <= 0.1,]$maxLRT.C)
minpC_elut <- min(AvsC_elut[AvsC_elut$qvalues.C > 0.1 & AvsC_elut$pvaluesFromqvalue_check.C < 0.01,]$maxLRT.C)

# for Pbra.A you need to do q value 0.2
minqA_pbra_0.2 <- min(AvsC_pbra[AvsC_pbra$qvalues.A <= 0.2,]$maxLRT.A)
minpA_pbra <- min(AvsC_pbra[AvsC_pbra$qvalues.A > 0.1 & AvsC_pbra$pvaluesFromqvalue_check.A < 0.01,]$maxLRT.A)
minqC_pbra <- min(AvsC_pbra[AvsC_pbra$qvalues.C <= 0.1,]$maxLRT.C)
minpC_pbra <- min(AvsC_pbra[AvsC_pbra$qvalues.C > 0.1 & AvsC_pbra$pvaluesFromqvalue_check.C < 0.01,]$maxLRT.C)

minqA_otters <- min(AvsC_otters[AvsC_otters$qvalues.A <= 0.1,]$maxLRT.A)
minpA_otters <- min(AvsC_otters[AvsC_otters$qvalues.A > 0.1 & AvsC_otters$pvaluesFromqvalue_check.A < 0.01,]$maxLRT.A)
minqC_otters <- min(AvsC_otters[AvsC_otters$qvalues.C <= 0.1,]$maxLRT.C)
minpC_otters <- min(AvsC_otters[AvsC_otters$qvalues.C > 0.1 & AvsC_otters$pvaluesFromqvalue_check.C < 0.01,]$maxLRT.C)

######## A vs. C, elut plot ######
plotAvsC_elut <- ggplot(AvsC_elut[AvsC_elut$Symbol.A!="N/A" & AvsC_elut$Symbol.C!="N/A",],aes(x=maxLRT.A,y=maxLRT.C,color=mergedLabel,shape=mergedLabel))+
  geom_point()+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_shape_manual(values=c(4,16,1))+
  scale_color_manual(values=c(artifactcol,realcol,unexaminedcol))+
  #geom_label_repel(data=subset(AvsC_elut,label=="putatively real"),aes(label=Symbol.C),size = 3,box.padding = 0.25,point.padding = 0.3,color="black")+
  geom_hline(yintercept = minqC_elut,color=filtercol1)+ # sig q value LRT cutoff for swamp 
  geom_vline(xintercept=minqA_elut,color=filtercol2)+ # sig q value LRT cutoff for gblocks
  geom_hline(yintercept = minpC_elut,color=filtercol1,linetype=2)+
  geom_vline(xintercept = minpA_elut,color=filtercol2,linetype=2)+
  geom_text(aes(30, minqC_elut+1),color=filtercol1, label = "q < 0.1", vjust = -1, size = 3.5)+
  geom_text(aes(30, minpC_elut+1),color=filtercol1, label = "p < 0.01", vjust = -1, size = 3.5)+
  geom_text(aes(minqA_elut+1,30),color=filtercol2, label = "q < 0.1",hjust=0,angle=90,size = 3.5)+
  geom_text(aes(minpA_elut+1,30),color=filtercol2, label = "p < 0.01", hjust=0,angle=90, size = 3.5)+
  ggtitle("Sea Otter Branch\nImpact of Stringent Swamp vs. Moderate Gblocks \nFiltered out genes with >1 BEB sites, Symbol = NA, Length < 120bp")+
  xlab("Stringent SWAMP maxLRT")+
  ylab("Moderate Gblocks maxLRT")

plotAvsC_elut

######## A vs. C, otters plot ######
plotAvsC_otters <- ggplot(AvsC_otters[AvsC_otters$Symbol.A!="N/A" & AvsC_otters$Symbol.C!="N/A",],aes(x=maxLRT.A,y=maxLRT.C,color=mergedLabel,shape=mergedLabel))+
  geom_point()+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_shape_manual(values=c(4,16,1))+
  scale_color_manual(values=c(artifactcol,realcol,unexaminedcol))+
  #geom_label_repel(data=subset(AvsC_otters,label=="putatively real"),aes(label=Symbol.C),size = 3,box.padding = 0.25,point.padding = 0.3,color="black")+
  geom_hline(yintercept = minqC_otters,color=filtercol1)+ # sig q value LRT cutoff for swamp 
  geom_vline(xintercept=minqA_otters,color=filtercol2)+ # sig q value LRT cutoff for gblocks
  geom_hline(yintercept = minpC_otters,color=filtercol1,linetype=2)+
  geom_vline(xintercept = minpA_otters,color=filtercol2,linetype=2)+
  geom_text(aes(25, minqC_otters+1),color=filtercol1, label = "q < 0.1", vjust = -1, size = 3.5)+
  geom_text(aes(25, minpC_otters-1),color=filtercol1, label = "p < 0.01", vjust = -1, size = 3.5)+
  geom_text(aes(minqA_otters+1,70),color=filtercol2, label = "q < 0.1",hjust=0,angle=90,size = 3.5)+
  geom_text(aes(minpA_otters+1,70),color=filtercol2, label = "p < 0.01", hjust=0,angle=90, size = 3.5)+
  ggtitle("Otters Branch\nImpact of Stringent Swamp vs. Moderate Gblocks \nFiltered out genes with >1 BEB sites, Symbol = NA, Length < 120bp")+
  xlab("Stringent SWAMP maxLRT")+
  ylab("Moderate Gblocks maxLRT")

plotAvsC_otters




#### A vs. C, pbra ####
plotAvsC_pbra <- ggplot(AvsC_pbra[AvsC_pbra$Symbol.A!="N/A" & AvsC_pbra$Symbol.C!="N/A",],aes(x=maxLRT.A,y=maxLRT.C,color=mergedLabel,shape=mergedLabel))+
  geom_point()+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_shape_manual(values=c(4,16,1))+
  scale_color_manual(values=c(artifactcol,realcol,unexaminedcol))+
  #geom_label_repel(data=subset(AvsC_pbra,label=="putatively real"),aes(label=Symbol.C),size = 3,box.padding = 0.25,point.padding = 0.3,color="black")+
  geom_hline(yintercept = minqC_pbra,color=filtercol1)+ # sig q value LRT cutoff for swamp 
  geom_vline(xintercept=minqA_pbra_0.2,color=filtercol2)+ # sig q value LRT cutoff for gblocks
  geom_hline(yintercept = minpC_pbra,color=filtercol1,linetype=2)+
  geom_vline(xintercept = minpA_pbra,color=filtercol2,linetype=2)+
  geom_text(aes(15, minqC_pbra+1),color=filtercol1, label = "q < 0.1", vjust = -1, size = 3.5)+
  geom_text(aes(15, minpC_pbra-1),color=filtercol1, label = "p < 0.01", vjust = -1, size = 3.5)+
  geom_text(aes(minqA_pbra_0.2+1,45),color=filtercol2, label = "q < 0.2",hjust=0,angle=90,size = 3.5)+
  geom_text(aes(minpA_pbra+1,45),color=filtercol2, label = "p < 0.01", hjust=0,angle=90, size = 3.5)+
  ggtitle("Giant Otter Branch\nImpact of Stringent Swamp vs. Moderate Gblocks \nFiltered out genes with >1 BEB sites, Symbol = NA, Length < 120bp")+
  xlab("Stringent SWAMP maxLRT")+
  ylab("Moderate Gblocks maxLRT")

plotAvsC_pbra

################### WRITE OUT GENE LISTS FOR TABLES ##########
########## ELUT GENE LIST TABLE #############
tableCols.1 = c("Symbol.A","TName","qvalues.A","pvaluesFromqvalue_check.A","qvalues.C","pvaluesFromqvalue_check.C")
tableCols.2 = c("Symbol","TName","qvalues","pvaluesFromqvalue_check")
# want all putatively real + p < 0.01 genes. Note that it's *OR* not *and* for the p < 0.01, because you want to have it be significant in EITHER, not just in both (mark those with both later on in the table)
AllElutGenesForTable.1 <- AvsC_elut[AvsC_elut$mergedLabel=="putatively real" & (AvsC_elut$pvaluesFromqvalue_check.A < 0.01 | AvsC_elut$pvaluesFromqvalue_check.C < 0.01),tableCols.1] # 61 genes
colnames(AllElutGenesForTable.1) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p","Gblocks_q","Gblocks_p")
# now need to add in those that didn't make it into others:
AllElutGenesForTable.2 <- elut_20170717_q[!(elut_20170717_q$Symbol %in% elut_20170922$Symbol) & elut_20170717_q$pvaluesFromqvalue_check < 0.01,tableCols.2] # 0 genes
colnames(AllElutGenesForTable.2) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p")
AllElutGenesForTable.2$Gblocks_q <- "x"
AllElutGenesForTable.2$Gblocks_p <- "x"

AllElutGenesForTable.3 <- elut_20170922_q[!(elut_20170922_q$Symbol %in% AllElutGenesForTable.1$Symbol) & elut_20170922_q$label=="putatively real" & elut_20170922_q$pvaluesFromqvalue_check < 0.01,tableCols.2] # 16 genes 
colnames(AllElutGenesForTable.3) <- c("Symbol","TranscriptID","Gblocks_q","Gblocks_p")
AllElutGenesForTable.3$Swamp_p <- "x"
AllElutGenesForTable.3$Swamp_q <- "x"

##### THESE ARE ALL ELUT GENES THAT ARE P < 0.01 IN EITHER 20170717 or 20170922
# and are putatively real in as well
AllElutGenesForTable.4 <- rbind(AllElutGenesForTable.1,AllElutGenesForTable.2,AllElutGenesForTable.3)
# check these are the same 
dim(AllElutGenesForTable.4)
length(unique(AllElutGenesForTable.4$Symbol))
write.table(AllElutGenesForTable.4,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/GeneTables/elut_20170717_20170922/elut.real.sigGenes.forMS.txt",row.names=F,quote=F,sep="\t")

########## PBRA GENE LIST TABLE #############
tableCols.1 = c("Symbol.A","TName","qvalues.A","pvaluesFromqvalue_check.A","qvalues.C","pvaluesFromqvalue_check.C")
tableCols.2 = c("Symbol","TName","qvalues","pvaluesFromqvalue_check")
# want all putatively real + p < 0.01 genes. Note that it's *OR* not *and* for the p < 0.01, because you want to have it be significant in EITHER, not just in both (mark those with both later on in the table)
AllPbraGenesForTable.1 <- AvsC_pbra[AvsC_pbra$mergedLabel=="putatively real" & (AvsC_pbra$pvaluesFromqvalue_check.A < 0.01 | AvsC_pbra$pvaluesFromqvalue_check.C < 0.01),tableCols.1] # 61 genes
colnames(AllPbraGenesForTable.1) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p","Gblocks_q","Gblocks_p")
dim(AllPbraGenesForTable.1) # 42 genes
# now need to add in those that didn't make it into others:
AllPbraGenesForTable.2 <- pbra_20170717_q[!(pbra_20170717_q$Symbol %in% pbra_20170922$Symbol) & pbra_20170717_q$pvaluesFromqvalue_check < 0.01,tableCols.2] # 0 genes
colnames(AllPbraGenesForTable.2) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p")
AllPbraGenesForTable.2$Gblocks_q <- "x"
AllPbraGenesForTable.2$Gblocks_p <- "x"
dim(AllPbraGenesForTable.2) # 6 genes (!)

AllPbraGenesForTable.3 <- pbra_20170922_q[!(pbra_20170922_q$Symbol %in% AllPbraGenesForTable.1$Symbol) & pbra_20170922_q$label=="putatively real" & pbra_20170922_q$pvaluesFromqvalue_check < 0.01,tableCols.2] # 10 genes 
colnames(AllPbraGenesForTable.3) <- c("Symbol","TranscriptID","Gblocks_q","Gblocks_p")
AllPbraGenesForTable.3$Swamp_p <- "x"
AllPbraGenesForTable.3$Swamp_q <- "x"
dim(AllPbraGenesForTable.3)
##### THESE ARE ALL PBRA GENES THAT ARE P < 0.01 IN EITHER 20170717 or 20170922
# and are putatively real in as well
AllPbraGenesForTable.4 <- rbind(AllPbraGenesForTable.1,AllPbraGenesForTable.2,AllPbraGenesForTable.3)
# check these are the same 
dim(AllPbraGenesForTable.4)
length(unique(AllPbraGenesForTable.4$Symbol))
write.table(AllPbraGenesForTable.4,"/Users/annabelbeichman/Documents/UCLA/Pbra/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/GeneTables/pbra_20170717_20170922/pbra.real.sigGenes.forMS.txt",row.names=F,quote=F,sep="\t")

########## OTTERS GENE LIST TABLE #############
tableCols.1 = c("Symbol.A","TName","qvalues.A","pvaluesFromqvalue_check.A","qvalues.C","pvaluesFromqvalue_check.C")
tableCols.2 = c("Symbol","TName","qvalues","pvaluesFromqvalue_check")
# want all putatively real + p < 0.01 genes. Note that it's *OR* not *and* for the p < 0.01, because you want to have it be significant in EITHER, not just in both (mark those with both later on in the table)
AllOttersGenesForTable.1 <- AvsC_otters[AvsC_otters$mergedLabel=="putatively real" & (AvsC_otters$pvaluesFromqvalue_check.A < 0.01 | AvsC_otters$pvaluesFromqvalue_check.C < 0.01),tableCols.1] #
colnames(AllOttersGenesForTable.1) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p","Gblocks_q","Gblocks_p")
dim(AllOttersGenesForTable.1)  #24
# now need to add in those that didn't make it into others:
AllOttersGenesForTable.2 <- otters_20170717_q[!(otters_20170717_q$Symbol %in% otters_20170922$Symbol) & otters_20170717_q$pvaluesFromqvalue_check < 0.01,tableCols.2] # 3 genes
colnames(AllOttersGenesForTable.2) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p")
AllOttersGenesForTable.2$Gblocks_q <- "x"
AllOttersGenesForTable.2$Gblocks_p <- "x"
dim(AllOttersGenesForTable.2) # 6 genes (!)

AllOttersGenesForTable.3 <- otters_20170922_q[!(otters_20170922_q$Symbol %in% AllOttersGenesForTable.1$Symbol) & otters_20170922_q$label=="putatively real" & otters_20170922_q$pvaluesFromqvalue_check < 0.01,tableCols.2] # 10 genes 
colnames(AllOttersGenesForTable.3) <- c("Symbol","TranscriptID","Gblocks_q","Gblocks_p")
AllOttersGenesForTable.3$Swamp_p <- "x"
AllOttersGenesForTable.3$Swamp_q <- "x"
dim(AllOttersGenesForTable.3) # 2 genes
##### THESE ARE ALL OTTERS GENES THAT ARE P < 0.01 IN EITHER 20170717 or 20170922
# and are putatively real in as well
AllOttersGenesForTable.4 <- rbind(AllOttersGenesForTable.1,AllOttersGenesForTable.2,AllOttersGenesForTable.3)
# check these are the same 
dim(AllOttersGenesForTable.4)
length(unique(AllOttersGenesForTable.4$Symbol))
write.table(AllOttersGenesForTable.4,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/GeneTables/otters_20170717_20170922/otters.real.sigGenes.forMS.txt",row.names=F,quote=F,sep="\t")

############# ALSO WRITE OUT ALL NON ARTIFACT + UNEXAMINED GENES for 20170717 for polysel ##########
## use the tables written out in step 8c instead (post VEP and NA filtering)
write.table(elut_20170717_q[elut_20170717_q$Symbol!="N/A" & elut_20170717_q$label!="artifact",],"output_codeml_20170717/inspection_elut/elut.putativelyReal.andUnexaminedGenes.20170717.forPolysel.txt",row.names=F,quote=F,sep="\t")
write.table(pbra_20170717_q[pbra_20170717_q$Symbol!="N/A" & pbra_20170717_q$label!="artifact",],"output_codeml_20170717/inspection_pbra/pbra.putativelyReal.andUnexaminedGenes.20170717.forPolysel.txt",row.names=F,quote=F,sep="\t")
write.table(otters_20170717_q[otters_20170717_q$Symbol!="N/A" & otters_20170717_q$label!="artifact",],"output_codeml_20170717/otters_inspection1/otters.putativelyReal.andUnexaminedGenes.20170717.forPolysel.txt",row.names=F,quote=F,sep="\t")
############################ VENN DIAGRAM CALCULATIONS ##################
###### elut: p < 0.01, not artifacts, intersections
## 20170922 vs 20170717 elut
gblocksSide_e1 <- elut_20170922_q[elut_20170922_q$Symbol!="N/A" & elut_20170922_q$uncorrPvalues.x < 0.01 & elut_20170922_q$label=="putatively real",]$Symbol
swampSide_e1 <- elut_20170717_q[elut_20170717_q$Symbol!="N/A" & elut_20170717_q$uncorrPvalues.x < 0.01 & elut_20170717_q$label!="artifact",]$Symbol

int_e1 <- intersect(gblocksSide_e1,swampSide_e1) # I didn't examine all the swamp genes if they were already ruled good in 20170922, I kept that designation
length(int_e1) # interior of venn diagram
write_clip(int_e1)
length(gblocksSide_e1)
length(swampSide_e1)
length(setdiff(gblocksSide_e1,swampSide_e1))
length(setdiff(swampSide_e1,gblocksSide_e1))

## 20170922 vs 20170717 pbra
gblocksSide_p1 <- pbra_20170922_q[pbra_20170922_q$Symbol!="N/A" & pbra_20170922_q$uncorrPvalues.x < 0.01 & pbra_20170922_q$label=="putatively real",]$Symbol
swampSide_p1 <- pbra_20170717_q[pbra_20170717_q$Symbol!="N/A" & pbra_20170717_q$uncorrPvalues.x < 0.01 & pbra_20170717_q$label!="artifact",]$Symbol

int_p1 <- intersect(gblocksSide_p1,swampSide_p1) # I didn't examine all the swamp genes if they were already ruled good in 20170922, I kept that designation
length(int_p1) # interior of venn diagram
write_clip(int_p1)
length(gblocksSide_p1)
length(swampSide_p1)
length(setdiff(gblocksSide_p1,swampSide_p1))
length(setdiff(swampSide_p1,gblocksSide_p1))

## 20170922 vs 20170717 otters
gblocksSide_o1 <- otters_20170922_q[otters_20170922_q$Symbol!="N/A" & otters_20170922_q$uncorrPvalues.x < 0.01 & otters_20170922_q$label=="putatively real",]$Symbol
swampSide_o1 <- otters_20170717_q[otters_20170717_q$Symbol!="N/A" & otters_20170717_q$uncorrPvalues.x < 0.01 & otters_20170717_q$label!="artifact",]$Symbol

int_o1 <- intersect(gblocksSide_o1,swampSide_o1) # I didn't examine all the swamp genes if they were already ruled good in 20170922, I kept that designation
length(int_o1) # interior of venn diagram
write_clip(int_o1)
length(gblocksSide_o1)
length(swampSide_o1)
length(setdiff(gblocksSide_o1,swampSide_o1))
length(setdiff(swampSide_o1,gblocksSide_o1))
p="\n")

################################ EXTRACT FOR VISUALIZATION ###################
cols_elut <- c("TName","group","cluster","sppCount","elut.bebCount","maxLRT","qvalues","pvaluesFromqvalue_check","Symbol")
cols_otters <- c("TName","group","cluster","sppCount","otters.bebCount","maxLRT","qvalues","pvaluesFromqvalue_check","Symbol")
cols_pbra <- c("TName","group","cluster","sppCount","pbra.bebCount","maxLRT","qvalues","pvaluesFromqvalue_check","Symbol")

###### 20170717: elut
elutInspect.20170717.pLT0.01 <- elut_20170717_q[elut_20170717_q$pvaluesFromqvalue_check < 0.01,]
write.table(elutInspect.20170717.pLT0.01[,cols],"output_codeml_20170717/inspection_elut/elutInspect.20170717.pLT0.01.txt",quote=F,row.names=F,col.names=T)
gatherSignificantResults_scriptMaker(sigResults = elutInspect.20170717.pLT0.01,branch="elut",date="20170717",note = "pLT0.01")
###### 20170922: elut
# inspect q > 0.1, p < 0.01 (what I think will be real)
elutInspect.20170922.qGT0.pLT0.01 <- elut_20170922_q[elut_20170922_q$qvalues >= 0.1 & elut_20170922_q$pvaluesFromqvalue_check < 0.01,]
# write out list for use in spreadsheet

write.table(elutInspect.20170922.qGT0.pLT0.01[,cols],"output_codeml_20170922/inspection_elut/elutInspect.20170922.qGT0.pLT0.01.txt",quote=F,row.names=F,col.names=T)
gatherSignificantResults_scriptMaker(sigResults = elutInspect.20170922.qGT0.pLT0.01,branch="elut",date="20170922",note = "qGT0.pLT0.01")

# inspect q < 0.1 # do spot check, don't have to check them all necessarily
elutInspect.20170922.qGT0.1 <- elut_20170922_q[elut_20170922_q$qvalues < 0.1,]
write.table(elutInspect.20170922.qGT0.1[,cols],"output_codeml_20170922/inspection_elut/elutInspect.20170922.qLT0.1.txt",quote=F,row.names=F,col.names=T)
gatherSignificantResults_scriptMaker(sigResults = elutInspect.20170922.qGT0.1,branch="elut",date="20170922",note = "qLT0.1")

## 20170922 otters extra that are new after convergence
newotters_20170922 <- otters_20170922_q[!(otters_20170922_q$Symbol %in% putativelyOkGenes_otters7$V1) & !(otters_20170922_q$Symbol %in% artifacts_otters7$V1) &  otters_20170922_q$pvaluesFromqvalue_check < 0.01,]

write.table(newotters_20170922[,cols_otters],"output_codeml_20170922/inspection_otters1/ottersInspect.20170922.newAfterConvergence.txt",quote=F,row.names=F,col.names=T)

gatherSignificantResults_scriptMaker(sigResults = newotters_20170922,branch="otters",date="20170922",note = "newAfterConvergence")

## 20170717 otters extra after convergence
newotters_20170717 <- otters_20170717_q[otters_20170717_q$pvaluesFromqvalue_check < 0.01 & !(otters_20170717_q$Symbol %in% artifacts_otters_20170717$Symbol) & !(otters_20170717_q$Symbol %in% putativelyOkGenes_otters_20170717$Symbol),]

write.table(newotters_20170717[,cols_otters],"output_codeml_20170717/otters_inspection1/ottersInspect.20170717.newAfterConvergence.txt",quote=F,row.names=F,col.names=T)

gatherSignificantResults_scriptMaker(sigResults = newotters_20170717,branch="otters",date="20170717",note = "newAfterConvergence")

## 20170921 otters extra after convergence
newotters_20170921 <- otters_20170921_q[otters_20170921_q$pvaluesFromqvalue_check < 0.01 & !(otters_20170921_q$Symbol %in% c(as.character(artifacts_otters6$Symbol),as.character(putativelyOkGenes_otters6$Symbol))),]

write.table(newotters_20170921[,cols_otters],"output_codeml_20170921/inspection_otters1/ottersInspect.20170921.newAfterConvergence.txt",quote=F,row.names=F,col.names=T)

gatherSignificantResults_scriptMaker(sigResults = newotters_20170921,branch="otters",date="20170921",note = "newAfterConvergence")


# 20170922 pbra with p LT 0.01 
pbraInspect.20170922.pLT0.01 <- pbra_20170922_q[pbra_20170922_q$pvaluesFromqvalue_check < 0.01,]
write.table(pbraInspect.20170922.pLT0.01[,cols_pbra],"output_codeml_20170922/inspection_pbra/pbraInspect.20170922.pLT0.01.txt",quote=F,row.names=F,col.names=T)
gatherSignificantResults_scriptMaker(sigResults = pbraInspect.20170922.pLT0.01,branch="pbra",date="20170922",note = "pLT0.01")




# 20170717 pbra with p LT 0.01 
pbraInspect.20170717.pLT0.01 <- pbra_20170717_q[pbra_20170717_q$pvaluesFromqvalue_check < 0.01,]
write.table(pbraInspect.20170717.pLT0.01[,cols_pbra],"output_codeml_20170717/inspection_pbra/pbraInspect.20170717.pLT0.01.txt",quote=F,row.names=F,col.names=T)
gatherSignificantResults_scriptMaker(sigResults = pbraInspect.20170717.pLT0.01,branch="pbra",date="20170717",note = "pLT0.01")
