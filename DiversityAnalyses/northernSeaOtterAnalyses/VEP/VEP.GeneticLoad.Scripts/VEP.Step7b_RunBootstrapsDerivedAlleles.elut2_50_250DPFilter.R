######## This is to run bootstraps on hoffman

### Files needed:
# 1. file with coding sites categories/# derived alleles (from step 7a)
# 2. file with average coding callable sites number 
args <- commandArgs(trailingOnly = TRUE) # have SGE_TASK_ID as 
# usage: Rscript $script ${SGE_TASK_ID}
bootID=args[1]
# do 10 at a time

# this was written out in step 6
# doing 100 groups of 10 bootstraps
nboot=10
spp="elut2"
date="20190215_nso_lib283only"
datadir=paste("/u/flashscratch/a/ab08028/otters/vcfs/vcf_",as.character(date),"_50_250_Filter/bootstrapDerivedAlleles/",sep="")
inputFile=paste(datadir,spp,".forBootstraps.StrippedDownDataset.AllCallableSites.derivedCountfromVEP.Annotation.50_250DPFilter.txt",sep="")
# get average coding callable from file: (need to put this file on Hoffman)
totalsCodingCallable <- read.table(paste(datadir,"/TotalCodingCallableSites.and.averageCodingCallableSiteselut1pbra.50_250DPFilter.elut2.lib283.txt",sep=""),header=F)
average_Coding_Callable <- totalsCodingCallable[totalsCodingCallable$V1=="average_Coding_Callable",]$V2
nsites=average_Coding_Callable # number of sites to draw out for each bootstrap (with replacement, use average callables sites in coding regions)
### Need to set up a dataframe of the right numbers of 
bootstrapDerivedAllelesCount <- function(inputFile,nboot,nsites){
  # set up empty results df:
  results <- data.frame(MissenseDerived=rep(NA,nboot),SynonymousDerived=rep(NA,nboot),Stop_GainedDerived=rep(NA,nboot),missenseHetCount=rep(NA,nboot),missenseHomCount=rep(NA,nboot),synHetCount=rep(NA,nboot),synHomGenotypeCount=rep(NA,nboot),sgHetGenotypeCount=rep(NA,nboot),sgHomCount=rep(NA,nboot))
  inputData <- read.table(inputFile,header=T)
  for(i in seq(1,nboot)){
    print(i)
    test <- inputData[sample(rownames(inputData),nsites,replace=TRUE),]
    # just do each one separately
    missenseCount <- sum(test[test$Annot=="missense",]$DerivedCount)
    synonymousCount <- sum(test[test$Annot=="synonymous",]$DerivedCount)
    stopgainedCount <- sum(test[test$Annot=="stopgained",]$DerivedCount)
    missenseHetCount <- sum(test$Annot=="missense" & test$DerivedCount==1)
    missenseHomCount <- sum(test$Annot=="missense" & test$DerivedCount==2)
    synHetCount <- sum(test$Annot=="synonymous" & test$DerivedCount==1)
    synHomCount <- sum(test$Annot=="synonymous" & test$DerivedCount==2)
    sgHetCount <- sum(test$Annot=="stopgained" & test$DerivedCount==1)
    sgHomCount  <- sum(test$Annot=="stopgained" & test$DerivedCount==2)
    results[i,] <- c(missenseCount,synonymousCount,stopgainedCount,missenseHetCount,missenseHomCount,synHetCount,synHomCount,sgHetCount,sgHomCount)
  }
  return(results)
}
# slow at first to read in file, then starts doing the bootstraps
results <- bootstrapDerivedAllelesCount(inputFile = inputFile,nboot=nboot,nsites=nsites)

# Then write out results
dir.create(paste(datadir,"/results/",sep=""))
write.table(results,paste(datadir,"/results/",spp,".group.",bootID,".",nboot,"bootstraps.results.txt",sep=""),row.names = F,col.names = T,quote=F,sep="\t")

