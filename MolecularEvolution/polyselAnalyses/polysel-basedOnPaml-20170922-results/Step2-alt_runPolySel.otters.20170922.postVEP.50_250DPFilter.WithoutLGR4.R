######################### Polysel run: 20170922 (3 reps, visual inspection pass ) ###############
## This is the latest script (20171219), where artifacts have been filtered out, and SO HAVE ~ 1400 genes flagged as having high impact VEP snps or indels based on alignment to ferret (in either elut or pbra)
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/CompGenomicsScripts/newPostCodemlAnalysis_DaubExcoffier_20170525/polysel")
# script and comments are following this guide: http://htmlpreview.github.io/?https://github.com/CMPG/polysel/blob/master/example_primates_homininae.html
# github: https://github.com/CMPG/polysel
require(qvalue)
require(Matrix)
require(igraph)
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#install.packages("Matrix")
#install.packages("igraph",dependencies = T)

branch="otters"
pamlDate="20170922" # gblocks, visual inspection pass

# set project variables
project.name=paste(branch,pamlDate,sep="_")
# format otters (project), paml.date, branch: (so could have paml.20170407/elut and /pbra too )
# 201701228: note that NAs are now included for characterized proteins that just lack a gene symbol (to be consistent with VEP)
datadir=paste("data/otters/paml.",pamlDate,"-traits-postVEP_50_250DPFilter-WithoutLGR4/",branch,sep="")
data.path=datadir
# format otters (project), paml.date, branch: (so could have paml.20170407/elut and /pbra too )
# set paths: make sure you've cloned github structure from Daub's github and its set up like hers
code.path<-"./R" # location of the R code for polysel
empfdr.path<-"./empfdr"
results.path<-"./results"

#The file polysel.R contains all the functions to run the pathway enrichment pipeline.
source(file.path(code.path,'polysel.R'))
### once you've set this stuff, the following code should be standard

######################## FILTER GENE SETS #################################
minsetsize<-10 # minimum # of genes in a set 
# changes it to remove genes that aren't in any sets obj.in.set=T
result<-ReadSetObjTables(in.path=data.path,
                         set.info.file="SetInfo.txt",
                         set.obj.file="SetObj.txt",
                         obj.info.file="ObjInfo.txt",
                         minsetsize=minsetsize,
                         obj.in.set=F,
                         merge.similar.sets=T)
# NOT excluding genes that aren't part of a set (good; too restrictive with so few sets)
set.info<-result$set.info
obj.info<-result$obj.info
set.obj<-result$set.obj
set.info.lnk<-result$set.info.lnk

cat("Number of sets: ", nrow(set.info), "\n", sep="")
cat("Number of genes: ", nrow(obj.info), "\n", sep="") 

print(set.info[1:5,],row.names=F, right=F)
print(obj.info[1:5,],row.names=F, right=F)
print(set.obj[1:5,],row.names=F, right=F)

########## CHECK FOR BIAS ######################
# They didn't end up having to correct for these things much. #
# can try some correlations: total otter gene length, used sequence length, tree length:
flds<-c("seqLength_mayInclGaps", "otters.treeLength","sppCount")
options(digits=2)
for (fld in flds){
  cat("Correlation between objStat and ", fld, ": ", sep="")
  ct<-(cor.test(obj.info$objStat,obj.info[[fld]]))
  cat("p-value: ", ct$p.value, " estimate: ", ct$estimate, 
      " [", ct$conf.int[1],
      ", ", ct$conf.int[2],"]\n",sep="")
  PlotStatField(obj.info, fld=fld, ylab=fld, xlab="objStat", 
                logaxis="y", show.bins=F)
}
### for now, not going to correct because correlations are very low (and similar to what Daub found and didn't correct for ) #######

#### CHECK SIZE DISTRUBTIONS ################
for (sz in c(10,50,250)){
  CheckStatDistribution(obj.info,setsize=sz,
                        n.rand = 100000, 
                        xlab="sum(objStat)")
}

# DAUB: After checking for biases and checking whether the null distribution is gaussian, and -if needed- correcting objStat for these biases, we create the data frame obj.stat which is derived from obj.info and consists of its fields objID, objStat and objBin. We save all the data objects so we can use them later on for testing.
######## SAVE DATA #############
obj.stat<-obj.info[,c("objID", "objStat", "objBin")]
save(set.info, obj.info, obj.stat, set.obj, set.info.lnk,
     file=file.path(data.path, "polysel_objects.RData"))

####### SET PARAMETERS ########
# DAUB:In our case, as we cannot approximate the SUMSTAT scores of random sets with a normal distribution, we have to create a null distribution by randomly sampling genes from the list of genes:
  
approx.null <- FALSE

#DAUB: Furthermore, we don’t assign the genes to bins: (only do this if you're correcting biases)

use.bins <- FALSE

#DAUB: If we have to create a null distribution (approx.null=FALSE), we use a sequential random sampling method to speed up the construction of null distributions. In short, for each tested gene set, we create a small null distribution and calculate a temporary p-value. For those sets that are below a certain p-value threshold, we repeatedly increase their corresponding null distribution until a maximum number of randomizations is reached.
seq.rnd.sampling <- TRUE

#DAUB: If we have to create a null distribution (approx.null=FALSE), we need to specify the number of random sets to create the null distribution to a sufficiently high value (> 400000 ). If you use sequential random sampling (seq.rnd.sampling = TRUE), you set the maximum number of randomizations. You can set to a lower value for nrand when you first want to test the function, but it should be larger than 1.


nrand=1000000


# DAUB: Usually we are looking for gene sets with extremely high SUMSTAT scores, so we compare their scores to the higher tail of the null distribution.

test <- "highertail"

#DAUB: Use lowertail or twosided if you are searching for low SUMSTAT scores or extreme scores respectively.

#DAUB: We use the function qvalue from the R package qvalue to calculate the q-value with the p-values of the tests.

#DAUB: We set here the parameter pi0.method of this function to bootstrap. See ?qvalue for more information about this function.

qvalue.method <- "bootstrap"

########## TEST WITHOUT PRUNING ############
# Test sets


result<-EnrichmentAnalysis(set.info, set.obj, obj.stat,
                           nrand=nrand, approx.null=approx.null, 
                           seq.rnd.sampling=seq.rnd.sampling,
                           use.bins=use.bins, test=test,
                           do.pruning=FALSE, minsetsize=minsetsize,
                           project.txt=project.name, do.emp.fdr=FALSE,
                           qvalue.method=qvalue.method
)

set.scores.prepruning <- result$set.scores.prepruning

#The ten highest scoring sets are (taking nrand = 1000):
  
print(set.scores.prepruning[1:10,],row.names=F, right=F)

save(set.scores.prepruning, 
     file = file.path(results.path,
                      paste(project.txt=project.name,
                            "_setscores_prepruning_", 
                            formatC(nrand,format="d"),".RData", sep="")))

write.table(set.scores.prepruning, quote=FALSE, sep="\t", row.names=FALSE, 
            file = file.path(results.path, paste(project.txt=project.name,                       "_setscores_prepruning_", formatC(nrand,format="d"),".txt",sep=""))) 



################ TEST WITH PRUNING #####################
result<-EnrichmentAnalysis(set.info, set.obj, obj.stat,
                           nrand=nrand, approx.null=approx.null, 
                           seq.rnd.sampling=seq.rnd.sampling,
                           use.bins=use.bins, test=test,
                           do.pruning=TRUE, minsetsize=minsetsize,
                           project.txt=project.name, do.emp.fdr=FALSE,
                           qvalue.method=qvalue.method
)

set.scores.postpruning <- result$set.scores.postpruning
print(set.scores.postpruning[1:10,],row.names=F, right=F)
####### go to the cluster now: Creating p-value null distribution on a CLUSTER ########## 
######### NOW MOVE OVER TO HOFFMAN (when it is back up) ##########
# on hoffman:
# do submit_step1.otters.paml.20170922-traits-postVEP_50_250DPFilter-WithLGR4.USETHISFORMANUSCRIPT.sh
# then submit_step2_polysel.runtests.20170922-traits-postVEP_50_250DPFilter-WithLGR4.USETHISFORMANUSCRIPT.sh


