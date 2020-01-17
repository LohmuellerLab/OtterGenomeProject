#### My version of Amber's olfaction cleaning script
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/elut1/")
######### get rid of duplicated scaffolds (once) ##########

blastsum <- read.table("step_1_ab_blastDb_blastResults/olfac_5sp_2seq_query.NOSPACES.sum",sep="/",header=T)
blastsumDD_raw <- blastsum

# First have to get rid of any hits to the duplicated scaffolds (this blasting was done before I deduped the genome)
# need to remove any queries where the Hit is in a dup scaff:
dupScaffs_99 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/DeDuplicateGenome_20170711/dup.scaffs.txt",header=F)

########## Want to pull out the duped hits 
blastsum_dup <- blastsum[(blastsum$Hit %in% dupScaffs_99$V1),]
dim(blastsum_dup) # 53982
#blastsumDD <- blastsum_dup
##### NOTE!!!! These lengths are in terms of Amino Acids but then when you do end - start it is nucleotides. So most should be > 810 NUCLEOTIDES (not AAs)
# try it w/out length filter
blastsumDD <- blastsum_dup[blastsum_dup$Hit_length >= 250,]

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

length(overlaps) 
head(overlaps)
# THIS IS STRANDED; so if they overlap but are on different strands they don't count
# get query numbers:
queryNums <- unique(queryHits(overlaps))

# "can use maxgap "intervals with a separation of maxgap or less and a minimum of minoverlap overlapping positions, allowing for maxgap, are considered to be overlapping. maxgap should be a scalar, non-negative, integer. minoverlap should be a scalar, positive integer." (Iranges manual)

# "In this case, and only this case, the drop.self and drop.redundant arguments are allowed. By default, the result will contain hits for each range against itself, and if there is a hit from A to B, there is also a hit for B to A. If drop.self is TRUE, all self matches are dropped. If drop.redundant is TRUE, only one of A->B and B->A is returned." (Iranges manual)

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
# these are the hits to duplicated scaffolds. So won't have any overlap with the genes I previously looked at. 
check <- findOverlaps(keep_prune_eval,maxgap=100,minoverlap=1,drop.self=F,drop.redundant=F) # 102 once length filter applied
# should only be self:self matches left! 

# write out
pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% keep_prune_eval$uniqueHitLabel,] # pull out the original entries 
############################# WRITE OUT .sum and .bed #####################################
write.table(pullOut,"cleanedByEvalue.LengthFilterOnly.Pruned.sum",row.names = F,quote=F,sep="/")

pullOut$bed_HitStart0 <- pullOut$Hit_Start - 1 # because bed is zero based (don't change end bc it is non inclusive)
# want to add 300 (100 AAs) to the start and end (doesn't matter what strand because adding in both directions)
pullOut$bed_HitStart0_minus300 <- pullOut$bed_HitStart0 - 300 # if goes negative, set to 0
pullOut$bed_HitStart0_minus300[pullOut$bed_HitStart0_minus300 < 0] <- 0
pullOut$bed_HitEnd_plus300 <- pullOut$Hit_END + 300 # some may be over limit of chromosome -- deal with that when it happens 
#### check to make sure you don't overhang:
# need other .fai ! 
scaffLengths <- read.table("sea_otter_23May2016_bS9RH.fasta.fai",header=F)
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
dim(bed) # 716
write.table(bed,"step_1_c_CleanedByEvalue_ABScript/DuplicateRemainder/elut1.DuplicateRemainderSeqs.cleanedByEvalue.250LengthFilter.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed",row.names = F,col.names = F,quote=F,sep="\t")


############################## Separate Clade Info ################
pullOut$clade <- sapply(strsplit(as.character(pullOut$Query), "_Clade"), "[", 2) # got this from: https://stackoverflow.com/questions/17215789/extract-a-substring-in-r-according-to-a-pattern
pullOut$ORspp <- sapply(strsplit(as.character(pullOut$Query), "_Clade"), "[", 1)

# add class I class I:
pullOut$class <- NA
pullOut[pullOut$clade=="_CLASS_I",]$class <- "Class I"
pullOut[pullOut$clade!="_CLASS_I",]$class <- "Class II"

write.table(pullOut,"step_1_c_CleanedByEvalue_ABScript/DuplicateRemainder/elut1.DuplicateRemainderSeq.cleanedByEvalue.250LengthFilter.Pruned.all.pullOut.info.txt",row.names = F,col.names = F,quote=F,sep="\t")
################### Plot #################################
require(ggplot2)
ggplot(pullOut[pullOut$Hit_length > 250,],aes(x=Hit,fill=class))+
  geom_histogram(stat="count")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 6))+
  ggtitle("Unique Blast Hits (removed overlaps)")

