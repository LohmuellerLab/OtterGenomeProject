#### My version of Amber's olfaction cleaning script
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/mfur/")
spp="mfur"
########### read back in the deduped (DD) results that you just wrote out (once) #########
blastsumDD_raw <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/mfur/step_1_ab_blastDb_blastResults/olfac_5sp_2seq_query.mfur.NOSPACES.sum",header=T,sep="/",quote="")
# do some filtering ahead of time to speed up processing:
# get rids of hits < 250bp in length
dim(blastsumDD_raw) # 411298 
##### NOTE!!!! These lengths are in terms of Amino Acids but then when you do end - start it is nucleotides. So most should be > 810 NUCLEOTIDES (not AAs)
### DON'T DO LENGTH CORRECTION FOR PGENE ANALYSIS 
blastsumDD <- blastsumDD_raw
dim(blastsumDD) # 411298 (318573 ONCE AAS < 250 REMOVED)

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
# for some reason, need that space between Query 1; Mamu whatever
head(grobj)
########## Find overlaps, maxgap=100,minoverlap=1,drop hits to itself, drop redundant hits #######
overlaps <- findOverlaps(grobj,maxgap=100,minoverlap=1,drop.self=F,drop.redundant=F) # don't need redundant hits 

length(overlaps) #

# get query numbers:
queryNums <- unique(queryHits(overlaps))
# 74627494 hits with strandedness 

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
# this can yield slightly diff output that original (non pgene analysis)
# if a <250bp hit has better Evalue than the original. Tends to only be 1-2 genes
# how to deal with: exclude closest (automate somehow? or just do by hand?)
keep_prune_eval <- cleanByEValue_prune(grobj,overlaps)

# check it! 
check <- findOverlaps(keep_prune_eval,maxgap=100,minoverlap=1,drop.self=F,drop.redundant=F)
# should only be self:self matches left! 
# THIS IS WORKING!!!
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
# need to do samtools faidx on genome:
scaffLengths <- read.table("../../20170707_MainOlfactionResults/mfur/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa.fai",header=F)
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
pullOut$bedname <- paste(spp,pullOut$Hit,pullOut$Hit_Start,pullOut$Hit_END,pullOut$Hit_strand,sep="_")
bed <- pullOut[,c("Hit","bed_HitStart0_minus300","bed_HitEnd_plus300","bedname","E.value","bedstrand")]
bed
dim(bed) # 716
write.table(bed,"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/mfur.cleanedByEvalue.noLengthFilter.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed",row.names = F,col.names = F,quote=F,sep="\t")
write.table(pullOut,"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/mfur.cleanedByEvalue.noLengthFilter.Pruned.all.pullOut.info.txt",row.names = F,col.names = T,quote=F,sep="\t")

############################## Separate Clade Info ################
pullOut$clade <- sapply(strsplit(as.character(pullOut$Query), "_Clade"), "[", 2) # got this from: https://stackoverflow.com/questions/17215789/extract-a-substring-in-r-according-to-a-pattern
pullOut$ORspp <- sapply(strsplit(as.character(pullOut$Query), "_Clade"), "[", 1)

# add class I class I:
pullOut$class <- NA
pullOut[pullOut$clade=="_CLASS_I",]$class <- "Class I"
pullOut[pullOut$clade!="_CLASS_I",]$class <- "Class II"

#### EXCLUDE FUNCTIONAL ORS  ####
pullOut <- read.table("step_1_c_CleanedByEvalue_ABScript_noLengthFilter/mfur.cleanedByEvalue.noLengthFilter.Pruned.all.pullOut.info.txt",sep="\t",header=T)
dim(pullOut)
head(pullOut)
head(pullOut)
functionalORs <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/mfur/step_9_FinalFunctionalList/classI.and.classII.final.ORs.passedInspection.txt",header=F) # this is the result of step 9 of functional analysis
dim(functionalORs)  # 816
length(unique(functionalORs$V1))
# so
functionalORs_names <- sapply(strsplit(as.character(functionalORs$V1),"-[0-9]+_aa",perl=T),"[",1)
length(unique(functionalORs_names))
# want to exclude any functional ORs
# want to exclude any functional ORs
nonFunctionalORs <- pullOut[!(pullOut$bedname %in% functionalORs_names),]
# if there are any functional ORs that aren't in pullOut (because of length stuff)
# exclude the one that beat them(? I guess?)
## CHECK
functionalORs_names[!(functionalORs_names %in% pullOut$bedname),] # missing 
# remove the ones that replaced them as they aren't really pgene, just an alternate chunk of the same
# View(pullOut)
# missing: "mfur_GL897276.1_595456_596208_-1" # how to deal with this missing dude. 
# find neighbor 
exclude <- "mfur_GL897276.1_595480_596211_-1"
dim(nonFunctionalORs) 
# 484 non functional
# make a bed file:
nonFunctionalORs <- nonFunctionalORs[!(nonFunctionalORs$bedname %in% exclude),]
dim(nonFunctionalORs) # 483 
nonFunctionalORs_bed <- nonFunctionalORs[,c("Hit","bed_HitStart0_minus300","bed_HitEnd_plus300","bedname","E.value","bedstrand")]
nonFunctionalORs_bed
write.table(nonFunctionalORs_bed,"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/nonFunctionalORs.mfur.bed",quote=F,row.names=F,col.names=F,sep="\t")


############################ WITHIN 30 of ends of SCAFF ################
# meant to be contigs not scaffs. Hm. 
# make a list of the ones that are close to ends of scaffs for now (contigs later?)

closeToStartOFScaffs <- nonFunctionalORs[nonFunctionalORs$Hit_Start < 30,]$bedname # close to start # 33 
closeToEndOFScaffs <- nonFunctionalORs[nonFunctionalORs$length - nonFunctionalORs$Hit_END <= 30,]$bedname # close to end of scaffold (39 regions)

write.table(c(as.character(closeToStartOFScaffs),as.character(closeToEndOFScaffs)),"step_1_c_CleanedByEvalue_ABScript_noLengthFilter/regionsCloseToStartEndofScaffold.mfur.txt",row.names=F,col.names=F,quote=F)

################### Plot #################################
require(ggplot2)
ggplot(pullOut,aes(x=Hit,fill=class))+
  geom_histogram(stat="count")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 6))+
  ggtitle("Unique Blast Hits (removed overlaps)")

