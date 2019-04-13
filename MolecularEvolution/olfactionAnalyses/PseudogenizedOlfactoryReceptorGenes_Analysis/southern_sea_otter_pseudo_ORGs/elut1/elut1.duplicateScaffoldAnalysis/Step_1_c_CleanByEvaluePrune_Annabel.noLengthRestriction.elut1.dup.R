#### My version of Amber's olfaction cleaning script
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut1/elut1.dup/")


########### read back in the deduped (DD) results that you just wrote out (once) #########
blastsum <- read.table("../../../20170707_MainOlfactionResults/elut1/step_1_ab_blastDb_blastResults/olfac_5sp_2seq_query.NOSPACES.sum",sep="/",header=T)
# do some filtering ahead of time to speed up processing:

##### NOTE!!!! These lengths are in terms of Amino Acids but then when you do end - start it is nucleotides. don't do the length filtering when looking for pseudogenes


dupScaffs_99 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/DeDuplicateGenome_20170711/dup.scaffs.txt",header=F)

blastsum_dup <- blastsum[(blastsum$Hit %in% dupScaffs_99$V1),]
dim(blastsum_dup) # 53982
blastsumDD <- blastsum_dup


dim(blastsumDD) # 53982 (217599 when removing < 250 AAs)

# want to add a unique index to each hit so that I can grab things out later.
blastsumDD$uniqueHitLabel <- seq(1:dim(blastsumDD)[1])

# 1. need to remove overlapping hits (within 100bp of each other) by choosing the one with the best e-value (lowest) or bit score (highest)
head(blastsumDD) 
# so what matters is Hit (scaffold ID), Hit_Start, Hit_End, Evalue and Bit Score of each hit.
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
require(GenomicRanges)

# seqnames are chromosomes/scaffolds
# ranges use IRanges to give starts then ends 
# strand you get from results
# score is bit score (or could be evalue)


grobj <- GRanges(seqnames=blastsumDD$Hit,ranges=IRanges(blastsumDD$Hit_Start,end=blastsumDD$Hit_END),strand=blastsumDD$Hit_strand,bit.score=blastsumDD$Bit.score,evalue=blastsumDD$E.value,query=blastsumDD$Query,hit.length=blastsumDD$Hit_length,uniqueHitLabel=blastsumDD$uniqueHitLabel)

head(grobj)
########## Find overlaps, maxgap=100,minoverlap=1,drop hits to itself, drop redundant hits #######
overlaps <- findOverlaps(grobj,maxgap=100,minoverlap=1,drop.self=F,drop.redundant=F) # don't need redundant hits 

length(overlaps) # withonly length filter: 74505009; without length filter: 149613979

queryNums <- unique(queryHits(overlaps))

# "can use maxgap "intervals with a separation of maxgap or less and a minimum of minoverlap overlapping positions, allowing for maxgap, are considered to be overlapping. maxgap should be a scalar, non-negative, integer. minoverlap should be a scalar, positive integer." (Iranges manual)

# "In this case, and only this case, the drop.self and drop.redundant arguments are allowed. By default, the result will contain hits for each range against itself, and if there is a hit from A to B, there is also a hit for B to A. If drop.self is TRUE, all self matches are dropped. If drop.redundant is TRUE, only one of A->B and B->A is returned." (Iranges manual)
# start with index 1 (group by index)
# for q in queryhits ...
# put in your genomic ranges object, and your overlaps result:

############################################ RUN PRUNING FUNCTION #######################
cleanByEValue_prune <- function(grobj,overlaps) {
  queryNums <- unique(queryHits(overlaps))
  keep <- GRanges()
  while(length(queryNums) > 0){
    q <- queryNums[1]
    print(q)
    subs_q <- subjectHits(overlaps[queryHits(overlaps)==q]) # don't need to add the query number in because I have the hits to itself kept into overlaps now     
    grobj_subs <- grobj[subs_q]
    keepThisHit <- which.min(grobj_subs$evalue)
    keep <- c(keep,grobj_subs[keepThisHit])
    queryNums <- queryNums[!(queryNums %in% subs_q)]
    print(paste("dropped:",subs_q))
    print(paste("left to go: ", length(queryNums)))
  }
  return(unique(keep))
}

keep_prune_eval <- cleanByEValue_prune(grobj,overlaps)

# check it! 
check <- findOverlaps(keep_prune_eval,maxgap=100,minoverlap=1,drop.self=F,drop.redundant=F) # note that all my final ORs are still in here, just have extras that would have failed length filter (checked it. 20171010)
check # 1114 up from 716 with length filter 
# should only be self:self matches left! 

# write out
pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% keep_prune_eval$uniqueHitLabel,] # pull out the original entries 
############################# WRITE OUT .sum and .bed #####################################
#write.table(pullOut,"cleanedByEvalue.LengthFilterOnly.Pruned.sum",row.names = F,quote=F,sep="/")
pullOut$bed_HitStart0 <- pullOut$Hit_Start - 1 # because bed is zero based (don't change end bc it is non inclusive)
# want to add 300 (100 AAs) to the start and end (doesn't matter what strand because adding in both directions)
pullOut$bed_HitStart0_minus300 <- pullOut$bed_HitStart0 - 300 # if goes negative, set to 0
pullOut$bed_HitStart0_minus300[pullOut$bed_HitStart0_minus300 < 0] <- 0
pullOut$bed_HitEnd_plus300 <- pullOut$Hit_END + 300 # some may be over limit of chromosome -- deal with that when it happens 
#### check to make sure you don't overhang:
#scaffLengths <- read.table("../../../20170707_MainOlfactionResults/elut1/sea_otter_23May2016_bS9RH.deduped.99.fasta.fai",header=F)
scaffLengths <- read.table("../sea_otter_23May2016_bS9RH.fasta.fai",header=F)

colnames(scaffLengths) <- c("Scaffold","length")
head(scaffLengths)
pullOut_lengths <- merge(pullOut,scaffLengths[,c("Scaffold","length")],by.x="Hit",by.y="Scaffold",all.x=T,all.y=F)
# check for overhangs:
pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]
# these are the bad ones that need to be fixed
pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]$bed_HitEnd_plus300 <- pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]$Hit_END

# check for overhangs again:
pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]
#none! good! 
pullOut <- pullOut_lengths


##### continue on:
# get strand:
pullOut$bedstrand <- "." # because bed needs + or - 
pullOut$bedstrand[pullOut$Hit_strand=="-1"] <- "-"
pullOut$bedstrand[pullOut$Hit_strand=="1"] <- "+"
# for the bedtools name want to have it be like Gang Li's script:
# $species\_$chr\_$start\_$end\_$string\ *** NOTE THAT BLAST IS ONE BASED *** so start end are 1 based
pullOut$bedname <- paste("elut",pullOut$Hit,pullOut$Hit_Start,pullOut$Hit_END,pullOut$Hit_strand,sep="_")
bed <- pullOut[,c("Hit","bed_HitStart0_minus300","bed_HitEnd_plus300","bedname","E.value","bedstrand")]
bed
dim(bed) # 204
#write.table(bed,"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/cleanedByEvalue.noLengthFilter.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed",row.names = F,col.names = F,quote=F,sep="\t")

######################## PSEUDOGENE SEARCH ####################
# UP til now this has been the same as step 1c from the functional analysis (except you got to keep < 250 AA hits)
# now want to read in your functional list:
#### EXCLUDE FUNCTIONAL ORS  #### 
functionalORs <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/elut1/elut1.duplicates/step_9_FinalFunctionalList/elut1.dup.classI.and.classII.final.ORs.passedInspection.txt",header=F) # this is the result of step 9 of functional analysis
dim(functionalORs)  # 61
# the names of these ORs have the extra info of length and "aa" at the end:
# elut_ScbS9RH_67127_851235_852047_1-348_aa
# want it to be
#elut_ScbS9RH_67127_851235_852047_1 to match your blast hits
# so
functionalORs_names <- sapply(strsplit(as.character(functionalORs$V1),"-3"),"[",1)
# want to exclude any functional ORs
nonFunctionalORs <- pullOut[!(pullOut$bedname %in% functionalORs_names),]
dim(nonFunctionalORs) 
# 649 non functional

# make a bed file:
nonFunctionalORs_bed <- nonFunctionalORs[,c("Hit","bed_HitStart0_minus300","bed_HitEnd_plus300","bedname","E.value","bedstrand")]
nonFunctionalORs_bed
write.table(nonFunctionalORs_bed,"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/nonFunctionalORs.bed",quote=F,row.names=F,col.names=F,sep="\t")

############################ WITHIN 30 of ends of SCAFF ################
# meant to be contigs not scaffs. Hm. 
# make a list of the ones that are close to ends of scaffs for now (contigs later?)

closeToStartOFScaffs <- nonFunctionalORs[nonFunctionalORs$Hit_Start < 30,]$bedname # close to start # 33 
closeToEndOFScaffs <- nonFunctionalORs[nonFunctionalORs$length - nonFunctionalORs$Hit_END <= 30,]$bedname # close to end of scaffold (39 regions)

write.table(c(closeToStartOFScaffs,closeToEndOFScaffs),"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/regionsCloseToStartEndofScaffold.txt",row.names=F,col.names=F,quote=F)





############################## Separate Clade Info ################
pullOut$clade <- sapply(strsplit(as.character(pullOut$Query), "_Clade"), "[", 2) # got this from: https://stackoverflow.com/questions/17215789/extract-a-substring-in-r-according-to-a-pattern
pullOut$ORspp <- sapply(strsplit(as.character(pullOut$Query), "_Clade"), "[", 1)

write.table(pullOut,"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/pullOut.cleanedByEvalue.noLengthFilter.Pruned.allHits.FuncAndNonFunc.txt",row.names = F,col.names = T,quote=F,sep="\t")


################### Plot #################################
require(ggplot2)
ggplot(pullOut,aes(x=Hit,fill=clade))+
  geom_histogram(stat="count")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 6))+
  ggtitle("Unique Blast Hits (removed overlaps)")

