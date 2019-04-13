setwd("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp")

################## STUFF YOU NEED TO RUN ##################
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

########### PACKAGES ###########
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
require(qvalue)
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
codemlGetPVals_allInfo_guidance_3Reps <- function(otters, elut, pbra,fdr=0.05){
  # for each spp. want to combine pvalues to do qvalue correction across all gene/branches
  allP <- rbind(otters[,c("branch","TName","uncorrPvalues")],elut[,c("branch","TName","uncorrPvalues")],pbra[,c("branch","TName","uncorrPvalues")])
  qobj <- qvalue(allP$uncorrPvalues, pi0.meth="bootstrap", fdr.level=fdr)
  allP$qvalues  <- qobj$qvalues
  allP$pvaluesFromqvalue_check <- qobj$pvalues
  allP$significant <- qobj$significant
  return(allP)
}
# this will give you info on each branch+TName
# then you can merge it with the others
########### test it: 20170922 (filtering method C)
otters_20170922 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170922,flag="otters",humanSymbols = hsaptt,dogSymbols = cfamtt)
elut_20170922 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170922,flag="elut",humanSymbols = hsaptt,dogSymbols = cfamtt)
pbra_20170922 <- codemlGetPVals_allInfo_guidance_3Reps(date=20170922,flag="pbra",humanSymbols = hsaptt,dogSymbols = cfamtt)

# do q value correction across all branches/genes:
qVals_20170922 <- codemlGetPVals_allInfo_guidance_3Reps(otters_20170922,elut_20170922,pbra_20170922,0.05)

# can then look at each by branch:
otters_20170922_q <- merge(otters_20170922,qVals_20170922[qVals_20170922$branch=="otters",],by="TName")
elut_20170922_q <- merge(elut_20170922,qVals_20170922[qVals_20170922$branch=="elut",],by="TName")
pbra_20170922_q <- merge(pbra_20170922,qVals_20170922[qVals_20170922$branch=="pbra",],by="TName")

