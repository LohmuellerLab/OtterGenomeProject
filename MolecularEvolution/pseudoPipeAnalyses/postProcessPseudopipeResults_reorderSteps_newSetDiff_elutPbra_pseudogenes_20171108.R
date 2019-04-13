############## think through sets
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/pseudoPipe_20170307")
require(clipr)
############################ GENES THAT ARE IN FERRET, NOT E,P ###############
# these are homologs that are in ferret and at least 4 other species, but are NOT in elut or pbra (and don't have paralogs in otehr spp.) These are pulled from my Portho results, so min identity is 25%
]# these are just the ferret protein IDs that met those criteria
hom_e <- read.table("pseudoPipeResults/ferretHomologs_notInElutPbraAnnotations/mfurHoms.min4Otherspp.notInElut.NoParalogs.txt",header=F) # these are in ferret but not elut (may or may not be in pbra)
head(hom_e)
colnames(hom_e) <- c("mfur","hsap","cfam")

hom_p <- read.table("pseudoPipeResults/ferretHomologs_notInElutPbraAnnotations/mfurHoms.min4Otherspp.notInPbra.NoParalogs.txt",header=F) # these are in ferret but not pbra (may or may not be in elut)
colnames(hom_p) <- c("mfur","hsap","cfam")

]

dim(hom_e) # 2669
dim(hom_p) # 2607

# A lot of these are probably missing due to incomplete annotation/assembly. However some may be pseudogenized! Let's see. 

###################################### MFUR: MFUR #######################################

# here is mfur peptides : mfur dna
mfur2mfur <- read.table("pseudoPipeResults/pseudoPipe_mfur_MfurPep_20170718/mfur_out_chunk_ppipe_MfurPep.allPgenes.20170718.txt",header=T)
head(mfur2mfur)
mfur2mfur$queryNo.1 <- gsub("\\.[0-9]+","",mfur2mfur$query) #



###################################### MFUR: ELUT #######################################
# here is mfur peptides: elut dna
# this has the headers of each run in it, but they are removed by comment character #
# top header I removed the # as part of the gather script
mfur2elut <- read.table("pseudoPipeResults/pseudoPipe_elut_MfurPep_20170717/enhydra_lutris_out_chunk_ppipe_MfurPep.allPgenes.20170717.txt",header=T)
### NOTE: was run on >15kb genome! With all gene models (including duplicate gene models on the duplicate scaffolds; so don't expect to see that many duplicated genes being miscalled as pseudogenes)
mfur2elut_dupScaffs99 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/DeDuplicateGenome_20170711/dup.scaffs.txt",header=F) # these are duplicate scaffods identified by dedupe with 99% identity. Use for mfur2elut. 
head(mfur2elut)
# Exclude those that hit to duplicated scaffolds:
# want to get rid of .1
mfur2elut$queryNo.1 <- gsub("\\.[0-9]+","",mfur2elut$query) #
# Exclude those that hit to duplicated scaffolds:
dim(mfur2elut) # 23227
############ REMOVED HITS TO DUPLICATED SCAFFS #########
mfur2elut <- mfur2elut[!(mfur2elut$chr %in% mfur2elut_dupScaffs99$V1),]
dim(mfur2elut) # 23149 (removed 78 hits)

###################################### MFUR: PBRA #######################################
mfur2pbra <- read.table("pseudoPipeResults/pseudoPipe_pbra_MfurPep_20170718/pteronura_brasiliensis_out_chunk_ppipe_MfurPep.allPgenes.20170718.txt",header=T)
head(mfur2pbra)
# want to get rid of .1
mfur2pbra$queryNo.1 <- gsub("\\.[0-9]+","",mfur2pbra$query) #


####################################### FILTER ##########################################################
# KEEP ALL TYPES TOGETHER FOR NOW 
##### 1. First exclude from pseudogenes pgenes that are found in elut/pbra as both a gene and a pseudogene (bc presumably still have function) ##### 
### these filters treat all 3 kinds of PGENE the same until the end
# intersect(mfur2elut$queryNo.1,hom_e$mfur) : these are genes that are pgenes in elut, are genes in ferret, and are not genes in elut (may or may not be genes in pbra)
elut_pgene1 <- mfur2elut[mfur2elut$queryNo.1 %in% intersect(mfur2elut$queryNo.1,hom_e$mfur),]
dim(elut_pgene1) # 3974
# these genes are pgene in elut, gene in ferret, not a gene in elut (Could be a gene in pbra)
pbra_pgene1 <- mfur2pbra[mfur2pbra$queryNo.1 %in% intersect(mfur2pbra$queryNo.1,hom_p$mfur),]
dim(pbra_pgene1) # 4950
# these genes are pgene in elut, gene in ferret, not a gene in elut (Could be a gene in elut)
# these two categories should both contain the ep genes (keep as 2 for now)

######  2. Now get rid of any pgenes that are also pgenes in mfur: ##### 
elut_pgene2 <- elut_pgene1[elut_pgene1$queryNo.1 %in% setdiff(elut_pgene1$queryNo.1,mfur2mfur$queryNo.1),]
# these are genes that are pgene in elut, gene in ferret, not a gene in elut and not a PGENE in mfur of any form 
dim(elut_pgene2) # 2188

pbra_pgene2 <- pbra_pgene1[pbra_pgene1$queryNo.1 %in% setdiff(pbra_pgene1$queryNo.1,mfur2mfur$queryNo.1),]
dim(pbra_pgene2) # 2993

##### 3 a. i. filter based on presence of disruptions: must not have 0 0 0 0 for ins del shift stop: must be a real disruption not just a fragment #####
# this is a new step as of Nov 2017
dim(elut_pgene2) # 2176
elut_pgene2 <- elut_pgene2[elut_pgene2$ins!=0 | elut_pgene2$del!=0 | elut_pgene2$shift!=0 | elut_pgene2$stop!=0, ]
dim(elut_pgene2) # 1538

dim(pbra_pgene2) # 2993
pbra_pgene2 <- pbra_pgene2[pbra_pgene2$ins!=0 | pbra_pgene2$del!=0 | pbra_pgene2$shift!=0 | pbra_pgene2$stop!=0, ]
dim(pbra_pgene2) # 1925
# (Now don't need to do this in step 4.)

######  3.b Categorize by whether these filtered pgenes are pgenes in elut only, pbra only, or in elut and pbra (should encompass all now) #####

elut_only_finalPgenes <- elut_pgene2[elut_pgene2$queryNo.1 %in% setdiff(elut_pgene2$queryNo.1,pbra_pgene2$queryNo.1),]
dim(elut_only_finalPgenes) # 572 now that >=0.75 restriciton (was 593 now that 0000s are removed) (was 799) (all three types)
length(unique(elut_only_finalPgenes$queryNo.1)) # 566  almost all are unique! 467 when unique (because there are hits that appear multiple times.)

pbra_only_finalPgenes <- pbra_pgene2[pbra_pgene2$queryNo.1 %in% setdiff(pbra_pgene2$queryNo.1,elut_pgene2$queryNo.1),]
dim(pbra_only_finalPgenes) # 338 once 0.75 removed (was 725 now that 0000s are removed (was 997))
length(unique(pbra_only_finalPgenes$queryNo.1)) # 322 unique (was 538 )

### NOW this must have disruptions in BOTH species (based on filter in step 3a). Good
# 20180102 --> took this out *AND* must be > 0.75 in both species 
ep_Pgenes_e <- elut_pgene2[elut_pgene2$queryNo.1 %in% intersect(elut_pgene2$queryNo.1,pbra_pgene2$queryNo.1),]
dim(ep_Pgenes_e) # 233 once 0.75 frac applied (was 945 now that any that is 0000 is removed.)  (was 1409) with multiple types
length(unique(ep_Pgenes_e$queryNo.1)) # 229 when unique (725) 

ep_Pgenes_p <- pbra_pgene2[pbra_pgene2$queryNo.1 %in% intersect(elut_pgene2$queryNo.1,pbra_pgene2$queryNo.1),]
dim(ep_Pgenes_p) # 233 
length(unique(ep_Pgenes_p$queryNo.1)) # 229 when unique 

#Interesting -- so are there more duplicates in pbra? That would make sense since it's more fragmented
length(unique(ep_Pgenes_p$queryNo.1)) # 928
# So these are the same length if duplicates removed; good. 
# Now I need to merge these:
ep_Pgenes_ep_FINAL <- merge(ep_Pgenes_e,ep_Pgenes_p,by="queryNo.1",suffixes=c(".elut",".pbra"))
dim(ep_Pgenes_ep_FINAL) # 243
# Number increases because there are some duplicates (e.g. a gene appears once in elut, but 3 x in Pbra. This then gets 3 entries in elut, that are identical in pbra.)
# this is the pbra info on the genes held in common
# then merge these with .elut and .pbra suffixes

# dim 3405 -- extras are from pgenes that are multiple types of pgene and so have several rows. If was unique in elut, but not in pbra, will fill in extra elut rows to account for that. 
# really just 928 unique sets that are in both elut and pbra. Some just appear as multiple entries. (likely to be annotation errors I suspect.)
############### venn diagram #########
#https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
#install.packages('VennDiagram')
library(VennDiagram)
grid.newpage()
## These are unique pseudogenes that have passed filters (are not both a gene and a psuedogene within the species, are a gene in ferret (+4 spp) and are not pseudogenized in ferret)
draw.pairwise.venn(area1 = length(unique(elut_pgene2$queryNo.1)), area2= length(unique(pbra_pgene2$queryNo.1)),cross.area = length(unique(ep_Pgenes_ep_FINAL$queryNo.1)),category = c("Elut","Pbra"),lty = rep("blank", 
2), fill = c("light blue", "light green"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
# ** note ** : some of these genes show up as multiple types of pseudogene, so there are many duplicates. 
##### 4. a. Get GConvert Gene Names and Descriptions ##### 
require(gProfileR)
# 4a. get gprofiler info
# AHA~ can do gconvert as part of gprofileR in R! woohoo!
elut_only_gconvert <- gconvert(unique(elut_only_finalPgenes$queryNo.1), organism = "mfuro", target = "ENSG", region_query = F,numeric_ns = "", mthreshold = Inf, filter_na = T, df = T)
dim(elut_only_gconvert)
pbra_only_gconvert <- gconvert(unique(pbra_only_finalPgenes$queryNo.1), organism = "mfuro", target = "ENSG", region_query = F,numeric_ns = "", mthreshold = Inf, filter_na = T, df = T)
# note: filter NA doesn't filter NAs from results, just gets rid of them if they are in input
dim(pbra_only_gconvert)
elut_and_pbra_gconvert <- gconvert(unique(ep_Pgenes_ep_FINAL$queryNo.1), organism = "mfuro", target = "ENSG", region_query = F,numeric_ns = "", mthreshold = Inf, filter_na = T, df = T)
dim(elut_and_pbra_gconvert)


#### # 4b. now want to merge this info with the pseudopipe info: ####
elut_only_finalPgenes_gconvert <- merge(elut_only_gconvert[,c("alias","target","name","description")],elut_only_finalPgenes,by.x="alias",by.y="queryNo.1")

pbra_only_finalPgenes_gconvert <- merge(pbra_only_gconvert[,c("alias","target","name","description")],pbra_only_finalPgenes,by.x="alias",by.y="queryNo.1")

ep_Pgenes_ep_FINAL_gconvert <- merge(elut_and_pbra_gconvert[,c("alias","target","name","description")],ep_Pgenes_ep_FINAL,by.x="alias",by.y="queryNo.1")

dim(elut_only_finalPgenes_gconvert) # 572


sum(unique(elut_only_finalPgenes_gconvert$name)!="N/A") # 493 (goes up slightly because move some frome shared to exclusive using filters) 401 (was 449 when all zeroes were allowed) unique genes that aren't N/A 
sum(unique(pbra_only_finalPgenes_gconvert$name)!="N/A") # 265 (464) (was 531) unique genes that aren't N/A
sum(unique(ep_Pgenes_ep_FINAL_gconvert$name)!="N/A") # 194 (642) (was 808) unique genes that aren't N/A 


#### 5. Separate by category -- optional. ####
## PSSD: 
sum(elut_only_finalPgenes_gconvert$type=="PSSD" & elut_only_finalPgenes_gconvert$name!="N/A") # 37
elut_only_finalPgenes_gconvert[elut_only_finalPgenes_gconvert$type=="PSSD" & elut_only_finalPgenes_gconvert$name!="N/A",]
sum(pbra_only_finalPgenes_gconvert$type=="PSSD" & pbra_only_finalPgenes_gconvert$name!="N/A") # 33
pbra_only_finalPgenes_gconvert[pbra_only_finalPgenes_gconvert$type=="PSSD" & pbra_only_finalPgenes_gconvert$name!="N/A",]
sum((ep_Pgenes_ep_FINAL_gconvert$type.elut=="PSSD" | ep_Pgenes_ep_FINAL_gconvert$type.pbra=="PSSD") & ep_Pgenes_ep_FINAL_gconvert$name!="N/A") # 56 where either elut or pbra has it as PSSD (other may have a different type)
ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$type.elut=="PSSD" & ep_Pgenes_ep_FINAL_gconvert$type.pbra=="PSSD" & ep_Pgenes_ep_FINAL_gconvert$name!="N/A",]


######### another way to categorize -- when both spp have them in common and it's the same type, examine:
ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$type.elut==ep_Pgenes_ep_FINAL_gconvert$type.pbra,]
length(unique(ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$type.elut==ep_Pgenes_ep_FINAL_gconvert$type.pbra,]$name)) # 761 are the same 

View(data.frame(unique(ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$type.elut==ep_Pgenes_ep_FINAL_gconvert$type.pbra,])))

############## VIEW ##########

View(pbra_only_finalPgenes_gconvert)
View(elut_only_finalPgenes_gconvert)
View(ep_Pgenes_ep_FINAL_gconvert)
View(mfur2mfur_gconvert)

############# 6. GPROFILER on frac > 0.75 ##############
### Need to fix version of Gprofiler. I am using the latest as of Nov 2017:
# http://biit.cs.ut.ee/gprofiler_archive3/r1732_e89_eg36/web/ 
# VERSION: r1732_e89_eg36 (November 9, 2017)

# older version more stable: try it
gProfileR::set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/r1730_e88_eg35/web/")
# elut-only
# pbra-only
# elut and pbra
# use homologs that aren't genes in elut as elut bg
# homologs that aren't genes in pbra as pbra bg

ferret="mfuro"
# ELUT:
GO_e_vsF <- gprofiler(as.vector(unique(as.character(elut_only_finalPgenes_gconvert[elut_only_finalPgenes_gconvert$frac >= 0.75,]$name))), organism = ferret,domain_size = "annotated",min_set_size = 2,max_set_size = 1000,min_isect_size = 3,src_filter = "HP",correction_method = "gSCS")
View(GO_e_vsF)
## add in shared pseudogenes because they could still interact biologically
# Ultimately decided to evaluate them separately to get at enrichment of traits through evolution, rather than present day(?)

GO_e_vsH <- (gprofiler(as.vector(unique(as.character(elut_only_finalPgenes_gconvert[elut_only_finalPgenes_gconvert$frac >= 0.75,]$name))), organism = "hsapiens",domain_size = "known",min_set_size = 2,max_set_size = 1000,min_isect_size = 3,src_filter = "HP",correction_method = "gSCS"))
dim(GO_e_vsH)
# PBRA:
GO_p_vsF <- (gprofiler(as.vector(unique(as.character(pbra_only_finalPgenes_gconvert[elut_only_finalPgenes_gconvert$frac >= 0.75,]$name))), organism = ferret,domain_size = "known",min_set_size = 2,max_set_size = 1000,min_isect_size = 3,src_filter = "HP",correction_method = "gSCS"))
View(GO_p_vsF)


GO_p_vsH <- (gprofiler(as.vector(unique(as.character(pbra_only_finalPgenes_gconvert[elut_only_finalPgenes_gconvert$frac >= 0.75,]$name))), organism = "hsapiens",domain_size = "known",min_set_size = 2,max_set_size = 1000,min_isect_size = 3,src_filter = "HP",correction_method = "gSCS"))
View(GO_p_vsH)


## Elut AND PBra: both fracs must be > 0.75

GO_ep_vsF <- (gprofiler(as.vector(unique(as.character(ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$frac.elut >= 0.75 & ep_Pgenes_ep_FINAL_gconvert$frac.pbra >= 0.75,]$name))), organism = ferret,domain_size = "known",min_set_size = 2,max_set_size = 1000,min_isect_size = 3,src_filter = "HP",correction_method = "gSCS",significant = T))
View(GO_ep_vsF)
# without those two genes you still get close p = 0.07, but no cigar.
# Still really interesting that those four genes are missing through -- even if not quite significant. Still indicates some hypoxia business.
GO_ep_vsH_no2genes <- (gprofiler(as.vector(unique(as.character(ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$frac.elut >= 0.75 & ep_Pgenes_ep_FINAL_gconvert$frac.pbra >= 0.75 & !(ep_Pgenes_ep_FINAL_gconvert$name %in% c("CASQ2","CACNA1C")),]$name))), organism = "hsapiens",domain_size = "known",min_set_size = 2,max_set_size = 1000,min_isect_size = 3,src_filter = "HP",correction_method = "gSCS",significant = T))
View(GO_ep_vsH_no2genes)  # Aha! Still significant against human background even once remove two genes that are dubious. How to systematically remove those genes though? Go by gene name in annotations?

# don't actually need >= 0.75 filter any more (did it in step 3a.i.)
GO_ep_vsH <- (gprofiler(as.vector(unique(as.character(ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$frac.elut >= 0.75 & ep_Pgenes_ep_FINAL_gconvert$frac.pbra >= 0.75,]$name))), organism = "hsapiens",domain_size = "known",min_set_size = 2,max_set_size = 1000,min_isect_size = 3,src_filter = "HP",correction_method = "gSCS"))
View(GO_ep_vsH)

######### WRITE OUT GO RESULTS ###########
write.table(GO_e_vsF,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_only/gProfilerResults/gProfiler.elut.FerretBG.frac0.75.txt",row.names = F,quote=F,sep="\t")
write.table(GO_e_vsF,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_only/gProfilerResults/gProfiler.elut.HumanBG.frac0.75.txt",row.names = F,quote=F,sep="\t")
write.table(GO_p_vsF,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/pbra_only/gProfilerResults/gProfiler.elut.FerretBG.frac0.75.txt",row.names = F,quote=F,sep="\t")
write.table(GO_p_vsF,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/pbra_only/gProfilerResults/gProfiler.elut.HumanBG.frac0.75.txt",row.names = F,quote=F,sep="\t")
# shared:
write.table(GO_ep_vsF,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_and_pbra/gProfilerResults/gProfiler.elutAndpbra.FerretBG.frac0.75.txt",row.names = F,quote=F,sep="\t")
write.table(GO_ep_vsH,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_and_pbra/gProfilerResults/gProfiler.elutAndpbra.HumanBG.frac0.75.txt",row.names = F,quote=F,sep="\t")


########### WRITE OUT SUPPLEMENTARY TABLES #####  
# Give all info of elut only, pbra only, ep_only genes. Each in own table with all the info. (Excel File for supplement)
# Make a nice table.
head(elut_only_finalPgenes_gconvert)
write.table(elut_only_finalPgenes_gconvert,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_only/elut_only_finalPgenes_gconvert_NonZeroDisruptions.noFracFilter.txt",quote=F,row.names=F,sep="\t")

#pbra:
write.table(pbra_only_finalPgenes_gconvert,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/pbra_only/pbra_only_finalPgenes_gconvert_NonZeroDisruptions.noFracFilter.txt",quote=F,row.names=F,sep="\t")
# shared: (wrong path! fixed it on 12/28/2017)
write.table(ep_Pgenes_ep_FINAL_gconvert,"pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_and_pbra/elut_and_pbra_shared_finalPgenes_gconvert_NonZeroDisruptions.noFracFilter.txt",quote=F,row.names=F,sep="\t")


######### Follow up on specific results ##############
## Get info on specific genes in GProfiler results:
ep_Pgenes_ep_FINAL_gconvert[as.character(ep_Pgenes_ep_FINAL_gconvert$name) %in% as.character(unlist(syncope_genes)),]

syncope_genes <- c("TECRL","CASQ2","THPO","HCRT","CACNA1C","P2RY11")
ep_Pgenes_ep_FINAL_gconvert[ep_Pgenes_ep_FINAL_gconvert$name %in% syncope_genes,]
## Retrieve sequences from genome.

