########### This script was used for final results and plots in the  ##########
## This script requires output from VEP analysis as well as codeml analysis (so you can exclude pseudogenized genes based on VEP

# The order goes: Step 8a >> VEP Steps 4a, 4b >> 8c [8b is an only-once step for lookup tables]. Final write out of 8c is what goes into polysel (has NA and VEP genes dealt with)
### Note that step 8a (previous version of this script is also important -- does the initial marking of things, but need this step too! Then it goes into Steb4b of VEP pipeline)
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
# I got the symbols from gene convert and manually added them using Google Sheets (not Excel! because of name conversions)

cfamtt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T)
hsaptt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T,quote="")
# dim(unique(hsaptt[,c("Transcript","Symbol")])) end up with 102915 unique combos, spread among more rows bc of alternate protein IDs etc

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
  # 20170925: adding 120bp length filter!
  ThreeReps_MAXlnls <- ThreeReps_MAXlnls[ThreeReps_MAXlnls$seqLength_mayInclGaps > 120,]
  # 20170925: dramatic move, but what If I get rid of beb >1 BEFORE qvalue?
  ThreeReps_MAXlnls <- ThreeReps_MAXlnls[ThreeReps_MAXlnls[,paste(flag,".bebCount",sep="")] <=1,]
  # convert to P value: df = 1
  ThreeReps_MAXlnls$uncorrPvalues <- (1-pchisq(ThreeReps_MAXlnls[,"maxLRT"],1))
  # want cluster info (names of other spp transcripts)
  # using CAT for this run because FNR was off in awk. In future use hsap.
  ThreeReps_MAXlnls_dogClusters <- merge(ThreeReps_MAXlnls,unique(dogSymbols[,c("cfam","Symbol")]),by.x="TName",by.y="cfam",all=F)
  ThreeReps_MAXlnls_humanClusters <- merge(ThreeReps_MAXlnls,unique(humanSymbols[,c("Transcript","Symbol")]),by.x="TName",by.y="Transcript",all=F)
  ThreeReps_MAXlnls_allClusters <- rbind(ThreeReps_MAXlnls_dogClusters,ThreeReps_MAXlnls_humanClusters)
  # this also gets all other spp and the gene symbol
  # want to get jbrowse 
  #return(codemldf_clusters)
  return(ThreeReps_MAXlnls_allClusters)
}
### try this! (need human and dog clusters)
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
##################### DEAL WITH N/A ISSUE ##############
# N/A genes didn't go through visual inspection; most are unchar proteins
# if any are not unchar proteins, inspect and leave in.
# elut, pbra, otters 20170717, 20170922, 20170921 >> symbol ==N/A, TName
# need to get rid of .1 after Tname for Uniprot
# Exclude any that are uncharacterized protein, keep the rest if they have any info
# look them up all together
elut_20170717_NA <- (lapply(strsplit(as.character(elut_20170717[elut_20170717$Symbol=="N/A",]$TName),"\\."),"[",1))
elut_20170922_NA <- (lapply(strsplit(as.character(elut_20170922[elut_20170922$Symbol=="N/A",]$TName),"\\."),"[",1))
pbra_20170717_NA <- (lapply(strsplit(as.character(pbra_20170717[pbra_20170717$Symbol=="N/A",]$TName),"\\."),"[",1))
pbra_20170922_NA <- (lapply(strsplit(as.character(pbra_20170922[pbra_20170922$Symbol=="N/A",]$TName),"\\."),"[",1))
otters_20170717_NA <- (lapply(strsplit(as.character(otters_20170717[otters_20170717$Symbol=="N/A",]$TName),"\\."),"[",1))
otters_20170921_NA <- (lapply(strsplit(as.character(otters_20170921[otters_20170921$Symbol=="N/A",]$TName),"\\."),"[",1))
otters_20170922_NA <- (lapply(strsplit(as.character(otters_20170922[otters_20170922$Symbol=="N/A",]$TName),"\\."),"[",1))

allNA <- unique(c(elut_20170717_NA,elut_20170922_NA,pbra_20170717_NA,pbra_20170922_NA,otters_20170922_NA,otters_20170921_NA,otters_20170717_NA)) # 352 genes

#write_clip(allNA) 

# Now go to Uniprot: http://www.uniprot.org/uploadlists/
# Copy paste in list
# do Ensembl Transcript >> Uniprot
# Download all results, read in:
uniprotNA <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/UniprotValidateNAs_allOutputs/uniprot_NA_results_elutPbraElut_TNameNASymbol_20171228",header=T,sep="\t")
dim(uniprotNA)
names(uniprotNA)
colnames(uniprotNA) <- c("yourlist","isomap","Entry", "Entry.name","Status","Protein.names" ,"Gene.names" ,"Organism" ,"Length","Annotation","Cross.reference..DisProt.")
# Any that are uncharacterized proteins, remove from dataset; keep rest in
UncharProteins <- uniprotNA[uniprotNA$Protein.names=="Uncharacterized protein",]
NotUnchar <- uniprotNA[uniprotNA$Protein.names!="Uncharacterized protein",]

### REMOVE SYMBOl=N/A and Uniprot is UNCHAR PROTEIN from datasets (uniprot) ####

######## 20170717

dim(elut_20170717_q)
elut_20170717_q$TName.No1 <- lapply(strsplit(as.character(elut_20170717_q$TName),"\\."),"[",1)
elut_20170717_q <- elut_20170717_q[!elut_20170717_q$TName.No1 %in% UncharProteins$yourlist,] ### removed uncharacterized proteins; left in symbol=NA that have a function in uniprot
dim(elut_20170717_q)

dim(pbra_20170717_q)
pbra_20170717_q$TName.No1 <- lapply(strsplit(as.character(pbra_20170717_q$TName),"\\."),"[",1)
pbra_20170717_q <- pbra_20170717_q[!pbra_20170717_q$TName.No1 %in% UncharProteins$yourlist,]
dim(pbra_20170717_q)

dim(otters_20170717_q)
otters_20170717_q$TName.No1 <- lapply(strsplit(as.character(otters_20170717_q$TName),"\\."),"[",1)
otters_20170717_q <- otters_20170717_q[!otters_20170717_q$TName.No1 %in% UncharProteins$yourlist,]
dim(otters_20170717_q)

######## 20170921
dim(elut_20170922_q)
elut_20170922_q$TName.No1 <- lapply(strsplit(as.character(elut_20170922_q$TName),"\\."),"[",1)
elut_20170922_q <- elut_20170922_q[!elut_20170922_q$TName.No1 %in% UncharProteins$yourlist,] ### removed uncharacterized proteins; left in symbol=NA that have a function in uniprot
dim(elut_20170922_q)

dim(pbra_20170922_q)
pbra_20170922_q$TName.No1 <- lapply(strsplit(as.character(pbra_20170922_q$TName),"\\."),"[",1)
pbra_20170922_q <- pbra_20170922_q[!pbra_20170922_q$TName.No1 %in% UncharProteins$yourlist,]
dim(pbra_20170922_q)

dim(otters_20170922_q)
otters_20170922_q$TName.No1 <- lapply(strsplit(as.character(otters_20170922_q$TName),"\\."),"[",1)
otters_20170922_q <- otters_20170922_q[!otters_20170922_q$TName.No1 %in% UncharProteins$yourlist,]
dim(otters_20170922_q)

######## 20170921
dim(otters_20170921_q)
otters_20170921_q$TName.No1 <- lapply(strsplit(as.character(otters_20170921_q$TName),"\\."),"[",1)
otters_20170921_q <- otters_20170921_q[!otters_20170921_q$TName.No1 %in% UncharProteins$yourlist,]
dim(otters_20170921_q)



####################### RESULTS OF VISUAL INSPECTION (includes NAs) ###########################
# as of 20171228 this includes NA symbols that are characterized but dont have a symbol, were inspected 
# for each branch and filter, inspect the genes w q < 0.1, and with p < 0.01 q > 0.1
####### vi: otters 20170921 #####
#artifacts_otters6 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170921/inspection_otters1/artifacts.20170921.txt",header=F)
artifacts_otters_20170921 <- read.table("output_codeml_20170921/inspection_otters1/artifacts_otters_20170921.txt",header=F)
colnames(artifacts_otters_20170921) <- "Symbol"
# didnt' do NA thing for 20170921, doesn't matter for now

#putativelyOkGenes_otters6 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170921/inspection_otters1/putativelyReal_20170921.txt",header=F)
putativelyOkGenes_otters_20170921 <- read.table("output_codeml_20170921/inspection_otters1/putatively_good_otters_20170921.txt",header=F)
colnames(putativelyOkGenes_otters_20170921) <- "Symbol"


# label the genes:
otters_20170921_q$label <- "unexamined"
otters_20170921_q[otters_20170921_q$Symbol %in% putativelyOkGenes_otters_20170921$Symbol,]$label <- "putatively real"
otters_20170921_q[otters_20170921_q$Symbol %in% artifacts_otters_20170921$Symbol,]$label <- "artifact"


################### vi: otters 20170922 ###############
# otters 20170922 (need to look at new genes after 3 reps, not too many though)
putativelyOkGenes_otters_20170922 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/putativelyReal_20170922.txt",header=F)
colnames(putativelyOkGenes_otters_20170922) <- "Symbol"
artifacts_otters_20170922 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/artifacts.20170922.txt",header=F)
colnames(artifacts_otters_20170922) <- "Symbol"

# note: otters 20170922 didn't have any NA artifact genes that needed to be removed  
# label the genes:
otters_20170922_q$label <- "unexamined"
otters_20170922_q[otters_20170922_q$Symbol %in% putativelyOkGenes_otters_20170922$Symbol,]$label <- "putatively real"
otters_20170922_q[otters_20170922_q$Symbol %in% artifacts_otters_20170922$Symbol,]$label <- "artifact"

################### vi: otters 20170717 ###############
putativelyOkGenes_otters_20170717 <- read.table("output_codeml_20170717/otters_inspection1/putativelyReal_20170717.txt",header=F)
colnames(putativelyOkGenes_otters_20170717) <- "Symbol"
artifacts_otters_20170717 <- read.table("output_codeml_20170717/otters_inspection1/artifacts.20170717.txt",header=F)
colnames(artifacts_otters_20170717) <- "Symbol"
artifacts_NA_otters_20170717 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/otters_inspection1/otters.20170717.NA.artifacts.txt") # these are the NA genes that aren't unchar proteins that were inspected and found to be artifacts (20171228)
# label the genes:
otters_20170717_q$label <- "unexamined"
otters_20170717_q[otters_20170717_q$Symbol %in% putativelyOkGenes_otters_20170717$Symbol,]$label <- "putatively real"
otters_20170717_q[otters_20170717_q$Symbol %in% artifacts_otters_20170717$Symbol,]$label <- "artifact"
otters_20170717_q[otters_20170717_q$TName %in% artifacts_NA_otters_20170717$V1,]$label <- "artifact"
################### vi: elut 20170922 ###############
# elut
putativelyOKGenes_elut_20170922 <- read.table("output_codeml_20170922/inspection_elut/elut.20170922.3reps.putativelyGOOD.txt",header=T)
artifacts_elut_20170922 <- read.table("output_codeml_20170922/inspection_elut/elut.20170922.3reps.ARTIFACT.txt",header=T)
artifacts_NA_elut_20170922 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_elut/elut.20170922.NA.Artifacts.txt")
# label the genes:
elut_20170922_q$label <- "unexamined"
elut_20170922_q[elut_20170922_q$Symbol %in% putativelyOKGenes_elut_20170922$Symbol,]$label <- "putatively real"
elut_20170922_q[elut_20170922_q$Symbol %in% artifacts_elut_20170922$Symbol,]$label <- "artifact"
elut_20170922_q[elut_20170922_q$TName %in% artifacts_NA_elut_20170922$V1,]$label <- "artifact"

################### vi: elut 20170717 ###############
# need to inspect the elut genes that are with p < 0.01 (both sig and non sig)
# for now only inspected the genes that weren't already inspected in 20170922's inspection (keeping real/artifact calls from 20170922) ***** <---- how did I deal with this??
putativelyOKGenes_elut_20170717 <- read.table("output_codeml_20170717/inspection_elut/elut.20170717.3reps.putativelyGOOD.txt",header=T) # note that these are only genes that aren't already inspected by 20170922
artifacts_elut_20170717 <- read.table("output_codeml_20170717/inspection_elut/elut.20170717.3reps.ARTIFACT.txt",header=T)
artifacts_NA_elut_20170717 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_elut/elut.20170717.NA.ARTIFACT.txt") # NA genes
# label the genes:
elut_20170717_q$label <- "unexamined"
elut_20170717_q[elut_20170717_q$Symbol %in% putativelyOKGenes_elut_20170717$Symbol,]$label <- "putatively real"
elut_20170717_q[elut_20170717_q$Symbol %in% artifacts_elut_20170717$Symbol,]$label <- "artifact"
elut_20170717_q[elut_20170717_q$TName %in% artifacts_NA_elut_20170717$V1,]$label <- "artifact"

#### 20180328: adding in the artifacts that I found when inspecting the 20170922 genes (didn't reinspect)
elut_20170717_q[elut_20170717_q$Symbol %in% artifacts_elut_20170922$Symbol,]$label <- "artifact"
elut_20170717_q[elut_20170717_q$Symbol %in% putativelyOKGenes_elut_20170922$Symbol,]$label <- "putatively real"

## 
################### vi: pbra 20170922 ###############
putativelyOkGenes_pbra_20170922 <- read.table("output_codeml_20170922/inspection_pbra/putatively_good_pbra_20170922.txt",header=F)
artifacts_pbra_20170922 <- read.table("output_codeml_20170922/inspection_pbra/artifacts_pbra_20170922.txt",header=F)
colnames(artifacts_pbra_20170922) <- "Symbol"
colnames(putativelyOkGenes_pbra_20170922) <- "Symbol"
artifacts_NA_pbra_20170922 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/pbra.20170922.NA.Artifacts.txt") # NA genes
putativelyOkGenes_NA_pbra_20170922 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/putatively_good_pbra_20170922.NA.txt")
# label the genes:
pbra_20170922_q$label <- "unexamined"
pbra_20170922_q[pbra_20170922_q$Symbol %in% putativelyOkGenes_pbra_20170922$Symbol,]$label <- "putatively real"
pbra_20170922_q[pbra_20170922_q$Symbol %in% artifacts_pbra_20170922$Symbol,]$label <- "artifact"
pbra_20170922_q[pbra_20170922_q$TName %in% artifacts_NA_pbra_20170922$V1,]$label <- "artifact"
pbra_20170922_q[pbra_20170922_q$TName=="ENSCAFT00000020689.3",]$label <- "artifact"
# **** note: failed to examine NA gene ENSCAFT00000020689 at first; did a quick look and looks like is caused by indels *** so added it to count; didn't update figure. # 
pbra_20170922_q[pbra_20170922_q$TName %in% putativelyOkGenes_NA_pbra_20170922$V1,]$label <- "putatively real" # only pbra 20170922 had an NA with p LT 0.01 that was putatively real

################### vi: pbra 20170717 ###############
putativelyOkGenes_pbra_20170717 <- read.table("output_codeml_20170717/inspection_pbra/putatively_good_pbra_20170717.txt",header=F)
artifacts_pbra_20170717 <- read.table("output_codeml_20170717/inspection_pbra/artifacts_pbra_20170717.txt",header=F)
artifacts_NA_pbra_20170717 <- read.table("~/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_pbra/pbra.20170717.NA.ARTIFACT.txt") # NA genes
colnames(artifacts_pbra_20170717) <- "Symbol"
colnames(putativelyOkGenes_pbra_20170717) <- "Symbol"
# label the genes:
pbra_20170717_q$label <- "unexamined"
pbra_20170717_q[pbra_20170717_q$Symbol %in% putativelyOkGenes_pbra_20170717$Symbol,]$label <- "putatively real"
pbra_20170717_q[pbra_20170717_q$Symbol %in% artifacts_pbra_20170717$Symbol,]$label <- "artifact"
pbra_20170717_q[pbra_20170717_q$TName %in% artifacts_NA_pbra_20170717$V1,]$label <- "artifact"
################################ RESULTS OF VEP ###########################
# these are the codeml clsuters flagged as having a SNP or indel that is HIGH impact in a protein DOMAIN in elut or pbra (in otters, if in elut OR pbra, got flagged )
#  (results of VEP Step 4a, 4b)
######### NEW: 20180207 I have updated this to be the results of VEP with 50/250% filter
# cds sequence only, CANONICAL, DOMAINS, HIGH IMPACT snps and indels!

########### vep: otters #########
# tehse are the results of VEP step4a and step4b!~ 
# old: 50/150 filter:
#otters_VEP_20170717 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/otters_inspection1/older_filterHighImpactVEP_50_150DPFilter/otters.VEPHighImpactGenes.gotRemoved.txt",sep="\t",header=T)
otters_VEP_20170717 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/otters_inspection1/otters.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",sep="\t",header=T)
# old: 
#otters_VEP_20170922 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/older_filterHighImpactVEP_50_150DPFilter/otters.VEPHighImpactGenes.gotRemoved.txt",sep="\t",header=T)
# new 50/250:
otters_VEP_20170922 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/otters.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",sep="\t",header=T)
# update labels:
otters_20170717_q$preVEPLabel <- otters_20170717_q$label
otters_20170717_q[otters_20170717_q$TName %in% otters_VEP_20170717$TName,]$label <- "VEP high impact"

otters_20170922_q$preVEPLabel <- otters_20170922_q$label
otters_20170922_q[otters_20170922_q$TName %in% otters_VEP_20170922$TName,]$label <- "VEP high impact"
########### vep: elut #########
# old: 
#elut_VEP_20170717 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_elut/older_filterHighImpactVEP_50_150DPFilter/elut.VEPHighImpactGenes.gotRemoved.txt",sep="\t",header=T)
# new: 
elut_VEP_20170717 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_elut/elut.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",sep="\t",header=T)
# old: 
#elut_VEP_20170922 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_elut/older_filterHighImpactVEP_50_150DPFilter/elut.VEPHighImpactGenes.gotRemoved.txt",sep="\t",header=T)
elut_VEP_20170922 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_elut/elut.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",sep="\t",header=T)
# update labels:
elut_20170717_q$preVEPLabel <- elut_20170717_q$label
elut_20170717_q[elut_20170717_q$TName %in% elut_VEP_20170717$TName,]$label <- "VEP high impact"

elut_20170922_q$preVEPLabel <- elut_20170922_q$label
elut_20170922_q[elut_20170922_q$TName %in% elut_VEP_20170922$TName,]$label <- "VEP high impact"

########### vep: pbra #########
# old:
#pbra_VEP_20170717 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_pbra/older_filterHighImpactVEP_50_150DPFilter/pbra.VEPHighImpactGenes.gotRemoved.txt",sep="\t",header=T)
pbra_VEP_20170717 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_pbra/pbra.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",sep="\t",header=T)
#old:
#pbra_VEP_20170922 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/older_filterHighImpactVEP_50_150DPFilter/pbra.VEPHighImpactGenes.gotRemoved.txt",sep="\t",header=T)
# new 
pbra_VEP_20170922 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/pbra.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",sep="\t",header=T)
# update labels:
pbra_20170717_q$preVEPLabel <- pbra_20170717_q$label
pbra_20170717_q[pbra_20170717_q$TName %in% pbra_VEP_20170717$TName,]$label <- "VEP high impact"

pbra_20170922_q$preVEPLabel <- pbra_20170922_q$label
pbra_20170922_q[pbra_20170922_q$TName %in% pbra_VEP_20170922$TName,]$label <- "VEP high impact"


########################### PLOTS FOR MANUSCRIPT: #######################

###################### COMPARE filter PLOTS ######################

########### OTTERS #######
AvsC_otters <- merge(otters_20170717_q,otters_20170922_q,by="TName",suffixes=c(".A",".C"))
# merge labels:
AvsC_otters$mergedLabel <- AvsC_otters$label.C
AvsC_otters[AvsC_otters$label.A=="artifact",]$mergedLabel <- "artifact"
AvsC_otters[AvsC_otters$label.A=="putatively real",]$mergedLabel <- "putatively real"
AvsC_otters[AvsC_otters$label.A=="VEP high impact",]$mergedLabel <- "VEP high impact" # this comes last because it trumps the visual inspection labels

########  ELUT ########
AvsC_elut <- merge(elut_20170717_q,elut_20170922_q,by="TName",suffixes=c(".A",".C"))
# merge labels:
AvsC_elut$mergedLabel <- AvsC_elut$label.C
AvsC_elut[AvsC_elut$label.A=="artifact",]$mergedLabel <- "artifact"
AvsC_elut[AvsC_elut$label.A=="putatively real",]$mergedLabel <- "putatively real"
AvsC_elut[AvsC_elut$label.A=="VEP high impact",]$mergedLabel <- "VEP high impact" # this comes last because it trumps the visual inspection labels
# PBRA
AvsC_pbra <- merge(pbra_20170717_q,pbra_20170922_q,by="TName",suffixes=c(".A",".C"))
# merge labels:
AvsC_pbra$mergedLabel <- AvsC_pbra$label.C
AvsC_pbra[AvsC_pbra$label.A=="artifact",]$mergedLabel <- "artifact"
AvsC_pbra[AvsC_pbra$label.A=="putatively real",]$mergedLabel <- "putatively real"
AvsC_pbra[AvsC_pbra$label.A=="VEP high impact",]$mergedLabel <- "VEP high impact" # this comes last because it trumps the visual inspection labels

# this is the minimum LRT value of genes in the category q < 0.1, or p < 0.01
minqA_elut <- min(AvsC_elut[AvsC_elut$qvalues.A <= 0.1,]$maxLRT.A)
minpA_elut <- min(AvsC_elut[AvsC_elut$qvalues.A > 0.1 & AvsC_elut$pvaluesFromqvalue_check.A < 0.01,]$maxLRT.A)
minqC_elut <- min(AvsC_elut[AvsC_elut$qvalues.C <= 0.1,]$maxLRT.C)
minpC_elut <- min(AvsC_elut[AvsC_elut$qvalues.C > 0.1 & AvsC_elut$pvaluesFromqvalue_check.C < 0.01,]$maxLRT.C)

######### for Pbra.A you need to do q value 0.2 #######
minqA_pbra_0.2 <- min(AvsC_pbra[AvsC_pbra$qvalues.A <= 0.2,]$maxLRT.A)
minpA_pbra <- min(AvsC_pbra[AvsC_pbra$qvalues.A > 0.1 & AvsC_pbra$pvaluesFromqvalue_check.A < 0.01,]$maxLRT.A)
minqC_pbra <- min(AvsC_pbra[AvsC_pbra$qvalues.C <= 0.1,]$maxLRT.C)
minpC_pbra <- min(AvsC_pbra[AvsC_pbra$qvalues.C > 0.1 & AvsC_pbra$pvaluesFromqvalue_check.C < 0.01,]$maxLRT.C)

minqA_otters <- min(AvsC_otters[AvsC_otters$qvalues.A <= 0.1,]$maxLRT.A)
minpA_otters <- min(AvsC_otters[AvsC_otters$qvalues.A > 0.1 & AvsC_otters$pvaluesFromqvalue_check.A < 0.01,]$maxLRT.A)
minqC_otters <- min(AvsC_otters[AvsC_otters$qvalues.C <= 0.1,]$maxLRT.C)
minpC_otters <- min(AvsC_otters[AvsC_otters$qvalues.C > 0.1 & AvsC_otters$pvaluesFromqvalue_check.C < 0.01,]$maxLRT.C)

######## A vs. C, elut plot ######
# used to be AvsC_elut[AvsC_elut$Symbol.A!="N/A" & AvsC_elut$Symbol.C!="N/A",]
# but now I have removed unchar proteins, so remaining NAs have some function, so want to keep them in an examine them (20171228)
plotAvsC_elut <- ggplot(AvsC_elut,aes(x=maxLRT.A,y=maxLRT.C,color=mergedLabel,shape=mergedLabel))+
  geom_point()+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_shape_manual(values=c(4,16,1,5))+
  scale_color_manual(values=c(artifactcol,realcol,unexaminedcol,artifactcol))+
  #geom_label_repel(data=subset(AvsC_elut,mergedLabel=="putatively real"),aes(label=Symbol.C),size = 3,box.padding = 0.25,point.padding = 0.3,color="black")+
  geom_hline(yintercept = minqC_elut,color=filtercol1)+ # sig q value LRT cutoff for swamp 
  geom_vline(xintercept=minqA_elut,color=filtercol2)+ # sig q value LRT cutoff for gblocks
  geom_hline(yintercept = minpC_elut,color=filtercol1,linetype=2)+
  geom_vline(xintercept = minpA_elut,color=filtercol2,linetype=2)+
  #geom_text(aes(30, minqC_elut+1),color=filtercol1, label = "q < 0.1", vjust = -1, size = 3.5)+
  #geom_text(aes(30, minpC_elut+1),color=filtercol1, label = "p < 0.01", vjust = -1, size = 3.5)+
  #geom_text(aes(minqA_elut+1,30),color=filtercol2, label = "q < 0.1",hjust=0,angle=90,size = 3.5)+
  #geom_text(aes(minpA_elut+1,30),color=filtercol2, label = "p < 0.01", hjust=0,angle=90, size = 3.5)+
  ggtitle("Sea Otter Branch\nImpact of Stringent Swamp vs. Moderate Gblocks \nFiltered out genes with >1 BEB sites, Symbol = NA, Length < 120bp")+
  xlab("SWAMP maxLRT")+
  ylab("Gblocks maxLRT")

plotAvsC_elut

######## A vs. C, otters plot ######
plotAvsC_otters <- ggplot(AvsC_otters,aes(x=maxLRT.A,y=maxLRT.C,color=mergedLabel,shape=mergedLabel))+
  geom_point()+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_shape_manual(values=c(4,16,1,5))+
  scale_color_manual(values=c(artifactcol,realcol,unexaminedcol,artifactcol))+
  #geom_label_repel(data=subset(AvsC_otters,label=="putatively real"),aes(label=Symbol.C),size = 3,box.padding = 0.25,point.padding = 0.3,color="black")+
  geom_hline(yintercept = minqC_otters,color=filtercol1)+ # sig q value LRT cutoff for swamp 
  geom_vline(xintercept=minqA_otters,color=filtercol2)+ # sig q value LRT cutoff for gblocks
  geom_hline(yintercept = minpC_otters,color=filtercol1,linetype=2)+
  geom_vline(xintercept = minpA_otters,color=filtercol2,linetype=2)+
  #geom_text(aes(30, minqC_otters+1),color=filtercol1, label = "q < 0.1", vjust = -1, size = 3.5)+
  #geom_text(aes(30, minpC_otters+1),color=filtercol1, label = "p < 0.01", vjust = -1, size = 3.5)+
  #geom_text(aes(minqA_otters+1,30),color=filtercol2, label = "q < 0.1",hjust=0,angle=90,size = 3.5)+
  #geom_text(aes(minpA_otters+1,30),color=filtercol2, label = "p < 0.01", hjust=0,angle=90, size = 3.5)+
  ggtitle("Otters Branch\nImpact of Stringent Swamp vs. Moderate Gblocks \nFiltered out genes with >1 BEB sites, Symbol = NA, Length < 120bp")+
  xlab("SWAMP maxLRT")+
  ylab("Gblocks maxLRT")

plotAvsC_otters




#### A vs. C, pbra ####
plotAvsC_pbra <- ggplot(AvsC_pbra,aes(x=maxLRT.A,y=maxLRT.C,color=mergedLabel,shape=mergedLabel))+
  geom_point()+
  theme_bw()+
  theme(legend.title=element_blank())+
  scale_shape_manual(values=c(4,16,1,5))+
  scale_color_manual(values=c(artifactcol,realcol,unexaminedcol,artifactcol))+
  #geom_label_repel(data=subset(AvsC_pbra,label=="putatively real"),aes(label=Symbol.C),size = 3,box.padding = 0.25,point.padding = 0.3,color="black")+
  geom_hline(yintercept = minqC_pbra,color=filtercol1)+ # sig q value LRT cutoff for swamp 
  geom_vline(xintercept=minqA_pbra_0.2,color=filtercol2)+ # sig q value LRT cutoff for gblocks
  geom_hline(yintercept = minpC_pbra,color=filtercol1,linetype=2)+
  geom_vline(xintercept = minpA_pbra,color=filtercol2,linetype=2)+
  #geom_text(aes(30, minqC_pbra+1),color=filtercol1, label = "q < 0.1", vjust = -1, size = 3.5)+
  #geom_text(aes(30, minpC_pbra+1),color=filtercol1, label = "p < 0.01", vjust = -1, size = 3.5)+
  #geom_text(aes(minqA_pbra+1,30),color=filtercol2, label = "q < 0.2",hjust=0,angle=90,size = 3.5)+
  #geom_text(aes(minpA_pbra+1,30),color=filtercol2, label = "p < 0.01", hjust=0,angle=90, size = 3.5)+
  ggtitle("Giant Otter Branch\nImpact of Stringent Swamp vs. Moderate Gblocks \nFiltered out genes with >1 BEB sites, Symbol = NA, Length < 120bp")+
  xlab("SWAMP maxLRT")+
  ylab("Gblocks maxLRT")

plotAvsC_pbra


############ save plots #########
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/compareFilters/postVEP/plotAvsC_elut.50_250DPFilter.pdf",plotAvsC_elut,device="pdf",width=7,height=5,units = "in")
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/compareFilters/postVEP/plotAvsC_pbra.50_250DPFilter.pdf",plotAvsC_pbra,device="pdf",width=7,height=5,units = "in")
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/compareFilters/postVEP/plotAvsC_otters.50_250DPFilter.pdf",plotAvsC_otters,device="pdf",width=7,height=5,units = "in")
################### WRITE OUT GENE LISTS FOR TABLES ##########

####################### TABLES FOR POLYSEL (USE THIS for polysel: 20171220) ##############
# want to write out non-artifact genes from 20170922 and 20170717 for polysel. only going to use 20170922
# != artifact, !=VEP-HIGH; symbol can be NA as long as unchar prots have been removed 

####### otters polysel tables ##########

# have to unlist TName.No1 column!
otters_20170717_q$TName.No1 <- unlist(otters_20170717_q$TName.No1)
writeOut_otters_20170717_q <- data.frame(otters_20170717_q[otters_20170717_q$label !="VEP high impact" & otters_20170717_q$label !="artifact",])
write.table(writeOut_otters_20170717_q,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/otters_inspection1/otters.fixedNAs.RemovedVEPHighImpactGenes.50_250DPFilter.putativelyReal.andUnexaminedGenes.20170717.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t") # instead, write this out in step 8c

# have to unlist TName.No1 column!
otters_20170922_q$TName.No1 <- unlist(otters_20170922_q$TName.No1)
writeOut_otters_20170922_q <- (otters_20170922_q[otters_20170922_q$label !="VEP high impact" & otters_20170922_q$label !="artifact",])
write.table(writeOut_otters_20170922_q,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/otters.fixedNAs.RemovedVEPHighImpactGenes.50_250DPFilter.putativelyReal.andUnexaminedGenes.20170922.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t") # instead, write this out in step 8c

####### pbra polysel tables ##########
pbra_20170717_q$TName.No1 <- unlist(pbra_20170717_q$TName.No1)
writeOut_pbra_20170717_q <- data.frame(pbra_20170717_q[pbra_20170717_q$label !="VEP high impact" & pbra_20170717_q$label !="artifact",])
write.table(writeOut_pbra_20170717_q,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_pbra/pbra.fixedNAs.RemovedVEPHighImpactGenes.50_250DPFilter.putativelyReal.andUnexaminedGenes.20170717.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t") # instead, write this out in step 8c

# have to unlist TName.No1 column!
pbra_20170922_q$TName.No1 <- unlist(pbra_20170922_q$TName.No1)
writeOut_pbra_20170922_q <- (pbra_20170922_q[pbra_20170922_q$label !="VEP high impact" & pbra_20170922_q$label !="artifact",])
write.table(writeOut_pbra_20170922_q,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/pbra.fixedNAs.RemovedVEPHighImpactGenes.50_250DPFilter.putativelyReal.andUnexaminedGenes.20170922.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t") # instead, write this out in step 8c

####### elut polysel tables ##########
elut_20170717_q$TName.No1 <- unlist(elut_20170717_q$TName.No1)
writeOut_elut_20170717_q <- data.frame(elut_20170717_q[elut_20170717_q$label !="VEP high impact" & elut_20170717_q$label !="artifact",])
write.table(writeOut_elut_20170717_q,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170717/inspection_elut/elut.fixedNAs.RemovedVEPHighImpactGenes.50_250DPFilter.putativelyReal.andUnexaminedGenes.20170717.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t") # instead, write this out in step 8c
# 20180301 fixed date in name ; right thing written out, name just wrong
# have to unlist TName.No1 column!
elut_20170922_q$TName.No1 <- unlist(elut_20170922_q$TName.No1)
writeOut_elut_20170922_q <- (elut_20170922_q[elut_20170922_q$label !="VEP high impact" & elut_20170922_q$label !="artifact",])
write.table(writeOut_elut_20170922_q,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_elut/elut.fixedNAs.RemovedVEPHighImpactGenes.50_250DPFilter.putativelyReal.andUnexaminedGenes.20170922.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t") # instead, write this out in step 8c


############## FOR MAIN TEXT / SUPPLEMENT (use this 20171228) ############
########## For supplement: ELUT GENE LIST TABLE #############
tableCols.1 = c("Symbol.A","TName","qvalues.A","pvaluesFromqvalue_check.A","qvalues.C","pvaluesFromqvalue_check.C")
tableCols.2 = c("Symbol","TName","qvalues","pvaluesFromqvalue_check")
# want all putatively real + p < 0.01 genes. Note that it's *OR* not *and* for the p < 0.01, because you want to have it be significant in EITHER, not just in both (mark those with both later on in the table)
AllElutGenesForTable.1 <- AvsC_elut[AvsC_elut$mergedLabel=="putatively real" & (AvsC_elut$pvaluesFromqvalue_check.A < 0.01 | AvsC_elut$pvaluesFromqvalue_check.C < 0.01),tableCols.1] # 61 genes
colnames(AllElutGenesForTable.1) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p","Gblocks_q","Gblocks_p")
# now need to add in those that didn't make it into others: *** THIS STEP IS BAD BECAUSE DOESN"T SPECIFY IT MUST BE NOT AN ARTIFACT. fixed it on 20180905; luckily not an issue for elut because there were none anyway.
AllElutGenesForTable.2 <- elut_20170717_q[!(elut_20170717_q$Symbol %in% elut_20170922$Symbol) & elut_20170717_q$pvaluesFromqvalue_check < 0.01 & elut_20170717_q$label=="putatively real",tableCols.2] # 0 genes
dim(AllElutGenesForTable.2)
colnames(AllElutGenesForTable.2) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p")
#AllElutGenesForTable.2$Gblocks_q <- "x" # don't need to do because no rows 
#AllElutGenesForTable.2$Gblocks_p <- "x"

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

write.table(AllElutGenesForTable.4,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/GeneTables/postVEP/elut_compare_20170717_20170922/elut.postVEP.real.sigGenes.forMS.50_250DPFilter.txt",row.names=F,quote=F,sep="\t") # checked this on 20180905 and it's okay for elut .

########## PBRA GENE LIST TABLE #############
tableCols.1 = c("Symbol.A","TName","qvalues.A","pvaluesFromqvalue_check.A","qvalues.C","pvaluesFromqvalue_check.C")
tableCols.2 = c("Symbol","TName","qvalues","pvaluesFromqvalue_check")
# want all putatively real + p < 0.01 genes. Note that it's *OR* not *and* for the p < 0.01, because you want to have it be significant in EITHER, not just in both (mark those with both later on in the table)
AllPbraGenesForTable.1 <- AvsC_pbra[AvsC_pbra$mergedLabel=="putatively real" & (AvsC_pbra$pvaluesFromqvalue_check.A < 0.01 | AvsC_pbra$pvaluesFromqvalue_check.C < 0.01),tableCols.1] # 61 genes
colnames(AllPbraGenesForTable.1) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p","Gblocks_q","Gblocks_p")
dim(AllPbraGenesForTable.1) # 42 genes
# now need to add in those that didn't make it into others: *** THIS STEP IS BAD BECAUSE DOESN"T SPECIFY IT MUST BE NOT AN ARTIFACT. fixed it on 20180905; luckily not an issue for elut because there were none anyway.
AllPbraGenesForTable.2 <- pbra_20170717_q[!(pbra_20170717_q$Symbol %in% pbra_20170922$Symbol) & pbra_20170717_q$pvaluesFromqvalue_check < 0.01 & pbra_20170717_q$label=="putatively real" ,tableCols.2] # 0 genes
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
write.table(AllPbraGenesForTable.4,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/GeneTables/postVEP/pbra_compare_20170717_20170922/pbra.postVEP.real.sigGenes.forMS.50_250DPFilter.txt",row.names=F,quote=F,sep="\t")

########## OTTERS GENE LIST TABLE #############
tableCols.1 = c("Symbol.A","TName","qvalues.A","pvaluesFromqvalue_check.A","qvalues.C","pvaluesFromqvalue_check.C")
tableCols.2 = c("Symbol","TName","qvalues","pvaluesFromqvalue_check")
# want all putatively real + p < 0.01 genes. Note that it's *OR* not *and* for the p < 0.01, because you want to have it be significant in EITHER, not just in both (mark those with both later on in the table)
AllOttersGenesForTable.1 <- AvsC_otters[AvsC_otters$mergedLabel=="putatively real" & (AvsC_otters$pvaluesFromqvalue_check.A < 0.01 | AvsC_otters$pvaluesFromqvalue_check.C < 0.01),tableCols.1] #
colnames(AllOttersGenesForTable.1) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p","Gblocks_q","Gblocks_p")
dim(AllOttersGenesForTable.1)  #24
# now need to add in those that didn't make it into others: *** THIS STEP IS BAD BECAUSE DOESN"T SPECIFY IT MUST BE NOT AN ARTIFACT. fixed it on 20180905; was not an issue for elut or pbra because there were no genes, but is an issue for otters: need to remove tmem8A and SLC35B3 from the manuscript; got it right in table S3.
AllOttersGenesForTable.2 <- otters_20170717_q[!(otters_20170717_q$Symbol %in% otters_20170922$Symbol) & otters_20170717_q$pvaluesFromqvalue_check < 0.01 & otters_20170717_q$label=="putatively real",tableCols.2] # 0 genes! 
colnames(AllOttersGenesForTable.2) <- c("Symbol","TranscriptID","Swamp_q","Swamp_p")
AllOttersGenesForTable.2$Gblocks_q <- "x"
AllOttersGenesForTable.2$Gblocks_p <- "x"
dim(AllOttersGenesForTable.2) # 20180905: 0 genes!!!! (previously 3 when doing it wrong) 

AllOttersGenesForTable.3 <- otters_20170922_q[!(otters_20170922_q$Symbol %in% AllOttersGenesForTable.1$Symbol) & otters_20170922_q$label=="putatively real" & otters_20170922_q$pvaluesFromqvalue_check < 0.01,tableCols.2] # 10 genes (2 genes?)
colnames(AllOttersGenesForTable.3) <- c("Symbol","TranscriptID","Gblocks_q","Gblocks_p")
AllOttersGenesForTable.3$Swamp_p <- "x"
AllOttersGenesForTable.3$Swamp_q <- "x"
dim(AllOttersGenesForTable.3) # 2 genes
##### THESE ARE ALL OTTERS GENES THAT ARE P < 0.01 IN EITHER 20170717 or 20170922
# and are putatively real in as well
AllOttersGenesForTable.4 <- rbind(AllOttersGenesForTable.1,AllOttersGenesForTable.2,AllOttersGenesForTable.3) # fixed this on 20180905
# check these are the same 
dim(AllOttersGenesForTable.4)
length(unique(AllOttersGenesForTable.4$Symbol))
write.table(AllOttersGenesForTable.4,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/newMSFigures_Tables/GeneTables/postVEP/otters_compare_20170717_20170922/otters.postVEP.real.sigGenes.forMS.50_250DPFilter.txt",row.names=F,quote=F,sep="\t") # overwrote old bad version of this on 20180905; will now update in manuscript



