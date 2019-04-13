#### Want to get pOrtho results
# and select out genes that have elut and/or pbra, and then some minimum number of other species. 
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/proteinOrtho_20170517_newspp")
portho <- read.table("15spp_POrtho_LongestIsosSelected_pruneTree_20170601/myproject.fixedHeader.poff",header=T,na.strings = c("*",","),stringsAsFactors = F)
head(portho)
head(portho,n=500)
### 1. replace any entries with commas with NA: 
portho_noCommas <- as.data.frame(sapply(portho,gsub,pattern=".*,.*",replacement=NA))
head(portho_noCommas,n=400)
### 2. Can't be NA in *both* elut and pbra (have to have at least one present)
portho_EP <- portho_noCommas[!is.na(portho_noCommas$elut.faa) | !is.na(portho_noCommas$pbra.faa),]
dim(portho_EP)
dim(portho)
dim(portho_noCommas)
### 3. Now want some limit on how many total species there should be. Going to set it to that there must be 8 (1 of which is elut or pbra; 1 of which is HSAP, gene must be known in human OR dog-- most characterized genomes, gives me a reliable ENSG ID, and if it isn't in human I find it's an uncharacterized protein). Need to recount now that I've gotten rid of commas. There are only 86 genes that are in at least 8 spp, but not in human or dog.
portho_EP$na_count <- apply(is.na(portho_EP), 1, sum)
head(portho_EP)
portho_EP$remainingSpp <- 15-portho_EP$na_count
head(portho_EP)
plot(hist(portho_EP$remainingSpp))
sum(portho_EP$remainingSpp >=8)
## must have human AND/OR dog: 
portho_EPHC <- portho_EP[!is.na(portho_EP$hsap.faa) | !is.na(portho_EP$cfam.faa),]
dim(portho_EPHC) # 15992
portho_EP_8 <- portho_EPHC[portho_EPHC$remainingSpp >=8,]
dim(portho_EP_8) # 15318 genes (# lost 86 genes that arent in cfam OR hsap)
sum(is.na(portho_EP_8$cfam.faa))
sum(is.na(portho_EP_8$hsap.faa)) # 593
sum(is.na(portho_EP$hsap.faa)) # 6038
sum(is.na(portho_EP_8$elut.faa))
sum(is.na(portho_EP_8$pbra.faa))
sum(is.na(portho_EP_8$mfur.faa))

# remember, using GUIDANCE2 you can drop more sequences that disrupt the alignments, this is just a preliminary filter
### 4. need to pick a name for each cluster: Going to use human and dog transcript IDs
head(portho_EP_8)
# first want to name them based on human, then on dog.
portho_EP_8$PName <- as.character(portho_EP_8$hsap.faa)

sum(is.na(portho_EP_8$PName)) # 715
portho_EP_8$PName[is.na(portho_EP_8$PName)] <- as.character(portho_EP_8$cfam.faa[is.na(portho_EP_8$PName)])


sppNum=15
#######  made lookup tables for each species: these give protein name: gene name : transcript ########

amelt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/amel.lookuptable.txt",stringsAsFactors = F)
colnames(amelt) <- c("Protein","Gene","Transcript")
# want an extra column that gives protein name without the .1 after each name (so it matches what I gave to portho)
amelt$ProteinNo.1 <- gsub("\\.[0-9]","",amelt$Protein)
head(amelt)

btaut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/btau.lookuptable.txt",stringsAsFactors = F)
colnames(btaut) <- c("Protein","Gene","Transcript")
btaut$ProteinNo.1 <- gsub("\\.[0-9]","",btaut$Protein)

cfamt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.lookuptable.txt",stringsAsFactors = F)
colnames(cfamt) <- c("Protein","Gene","Transcript")
cfamt$ProteinNo.1 <- gsub("\\.[0-9]","",cfamt$Protein)

ecabt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/ecab.lookuptable.txt",stringsAsFactors = F)
colnames(ecabt) <- c("Protein","Gene","Transcript")
ecabt$ProteinNo.1 <- gsub("\\.[0-9]","",ecabt$Protein)

fcatt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/fcat.lookuptable.txt",stringsAsFactors = F)
colnames(fcatt) <- c("Protein","Gene","Transcript")
fcatt$ProteinNo.1 <- gsub("\\.[0-9]","",fcatt$Protein)

hsapt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.lookuptable.txt",stringsAsFactors = F)
colnames(hsapt) <- c("Protein","Gene","Transcript")
hsapt$ProteinNo.1 <- gsub("\\.[0-9]","",hsapt$Protein)

mfurt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/mfur.lookuptable.txt",stringsAsFactors = F)
colnames(mfurt) <- c("Protein","Gene","Transcript")
mfurt$ProteinNo.1 <- gsub("\\.[0-9]","",mfurt$Protein)

ttrut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/proteinOrtho_20170517_newspp/longestIsoform_ttru_pvam_mluc/ttru/ttru.lookuptable.txt",stringsAsFactors = F)
colnames(ttrut) <- c("Protein","Gene","Transcript")
ttrut$ProteinNo.1 <- gsub("\\.[0-9]","",ttrut$Protein)

pvamt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/proteinOrtho_20170517_newspp/longestIsoform_ttru_pvam_mluc/pvam/pvam.lookuptable.txt",stringsAsFactors = F)
colnames(pvamt) <- c("Protein","Gene","Transcript")
pvamt$ProteinNo.1 <- gsub("\\.[0-9]","",pvamt$Protein)


mluct <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/proteinOrtho_20170517_newspp/longestIsoform_ttru_pvam_mluc/mluc/mluc.lookuptable.txt",stringsAsFactors = F)
colnames(mluct) <- c("Protein","Gene","Transcript")
mluct$ProteinNo.1 <- gsub("\\.[0-9]","",mluct$Protein)
# dont need umart, orost, lwedt lookuptables here, becuase protein and transcript have same name

## can combine each spp look up table into one big lookup table:
lookUpAll <- rbind(amelt,btaut,cfamt,ecabt,fcatt,hsapt,mfurt,ttrut,pvamt,mluct)
dim(lookUpAll)
# add to lookup table: if protein is NA, put in NA
head(lookUpAll)
lookUpAll2 <- rbind(lookUpAll,c(NA,NA,NA,NA)) # need to add this so it knows to put an NA where the TName is
########### protein:transcript ###############

# this assumes lookup table is as above: protein, gene, transcript, proteinID with out .1 after it
# and skips the first four columns of all the one2one POrthoResults 
# requires dplyr
library(dplyr) ### note: Dplyr and plyr are incompatible
GetTranscriptIDs <- function(spp,lookuptable,one2onePOrthos){
  print(spp)
  # merge spp transcript IDs in:
  merged <- merge(one2onePOrthos,lookuptable[,3:4],by.x=paste(spp,".faa",sep=""),by.y="ProteinNo.1")
  return(merged)
}
# dont need a sea otter  or giant otter look up table; just need to use otter names from annotation;same for lwed, umar, oros
animals <- c("amel","btau","cfam","ecab","fcat","hsap","mfur","ttru","pvam","mluc")
t1 <- portho_EP_8[,4:length(portho_EP_8)]
library(data.table) # requires this library for setnames
for(animal in animals){
  t1 <- GetTranscriptIDs(animal,lookUpAll2,t1)
  setnames(t1,old="Transcript",new=paste(animal,"_T",sep=""))
}

head(t1)
dim(t1) # 15318 19 , should be 15318
# want to keep the out elut.faa (sea otter) and pbra.faa (giant otter) umar, oros, lwed.faa columns, and all the _T cols:
# want to add a prefix to every lwed, oros, umar entry so that they can be distinguished later. # note this also alters the NA entries
t1$oros.faa <- paste("OROS",t1$oros.faa,sep="_")
head(t1)
t1$lwed.faa <- paste("LWED",t1$lwed.faa,sep="_")
head(t1)
t1$umar.faa <- paste("UMAR",t1$umar.faa,sep="_")
head(t1)
# and discard the protein IDs for the other species. 
tt <- t1[,(sppNum-4):length(t1)] ## need to change the -4 depending on # of species!
# get rid of na_count, remaining spp and PName columns:
tt <- subset(tt, select=-c(na_count,remainingSpp,PName))
# need to give umar, oros and lwed a header 
head(tt) ### make sure this has all the species! 
ttdf <- t(as.data.frame(tt,stringsAsFactors = F))
#check it: first 10 clusters
ttdf[,1:10]
dim(ttdf) # 12 9309
colnames(ttdf) <- c(seq(1:dim(ttdf)[2]))
ttdf[,1:10]
### this ttdf contains each cluster in 'wide' format
# so each column is a different gene cluster
# and each row is a different species

# write out table of ortholog clusters passing filters:
write.table(ttdf,"15spp_POrtho_LongestIsosSelected_pruneTree_20170601/one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt",col.names = F, quote = F,row.names=F)

