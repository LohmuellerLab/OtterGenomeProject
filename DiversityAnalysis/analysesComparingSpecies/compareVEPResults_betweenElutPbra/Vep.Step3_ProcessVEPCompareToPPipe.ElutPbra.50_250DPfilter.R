require(gProfileR)
require(ggplot2)
require(clipr)
#source("https://bioconductor.org/biocLite.R")
#biocLite("GOsummaries")
#require(GOsummaries)

############ ELUT ##################
####### 1. read in VEP results #######

#### look at HOMOZYGOUS alt only: 
elut_vep_indels_high_domain <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/01_Elut_CA_Gidget.raw_variants.cdsSequence.AllScaffsConcat.20171006.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
# these are hom. alt. INDELS
colnames(elut_vep_indels_high_domain) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(elut_vep_indels_high_domain) # 5167 (previously was 4992 (~80 are stopgained))
length(unique(elut_vep_indels_high_domain$Gene)) #2373 unique genes (gained 90) (was 2283 unique genes (not all were in paml clusters, that is why numbers are different)) 

# new:
elut_vep_snps_high_domain <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/01_Elut_CA_Gidget.raw_variants.cdsSequence.AllScaffsConcat.20171006.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
colnames(elut_vep_snps_high_domain) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(elut_vep_snps_high_domain) # 711 (was 679)
length(unique(elut_vep_snps_high_domain$Gene)) # 642 distinct genes (was 609)
############ get all high impact veps #########
elut_vep_snps_high_domain$category <- "SNP"
elut_vep_indels_high_domain$category <- "Indel"
elut_allHighDomain <- rbind(elut_vep_indels_high_domain,elut_vep_snps_high_domain)
dim(elut_allHighDomain)
############## get gene IDs ###############

elut_allHighDomain_nm <- gconvert(query=elut_allHighDomain$Gene,organism="mfuro")
dim(elut_allHighDomain_nm)


############ PBRA ##################
####### 1. read in PBRA results #######


pbra_vep_indels_high_domain <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/PteBra_1.raw_variants.cdsSequence.AllScaffsConcat.20171206.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
colnames(pbra_vep_indels_high_domain) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(pbra_vep_indels_high_domain) # 4544
length(unique(pbra_vep_indels_high_domain$Gene)) # 2180

pbra_vep_snps_high_domain <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/PteBra_1.raw_variants.cdsSequence.AllScaffsConcat.20171206.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
colnames(pbra_vep_snps_high_domain) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(pbra_vep_snps_high_domain) # 705
length(unique(pbra_vep_snps_high_domain$Gene)) # 623 distinct genes
############ get all high impact veps #########
pbra_vep_snps_high_domain$category <- "SNP"
pbra_vep_indels_high_domain$category <- "Indel"
pbra_allHighDomain <- rbind(pbra_vep_indels_high_domain,pbra_vep_snps_high_domain)
dim(pbra_allHighDomain)
############## get gene IDs ###############
pbra_allHighDomain_nm <- gconvert(query=pbra_allHighDomain$Gene,organism="mfuro")
head(pbra_allHighDomain_nm)
# figure out which NA symbol genes are uncharacterized:


######################## LOOK UP GENES IN UNIPROT, discard  "uncharacterized protein" #######

write_clip(pbra_allHighDomain_nm[pbra_allHighDomain_nm$name=="N/A",]$alias)
# then go to uniprot:
# http://www.uniprot.org/uniprot/ search and save results to
# 50/150 results: pbraUniprot <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP_50_150_ResultsScripts/Steps3-4-CompareElutPbra/Uniprot_VEP_Results/uniprotVEP_onlineSearch_20171228/pbra.VEP_NAGenes.UniprotResults.20171228",header=T,sep="\t")
# 50/250 results:
pbraUniprot <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/uniprotVEP_onlineSearch_20180207/pbra.VEP_NAGenes.UniprotResults.50_250DP.20180207",header=T,sep="\t")
head(pbraUniprot)
colnames(pbraUniprot) <- c("yourlist","isomap","Entry", "Entry.name","Status","Protein.names" ,"Gene.names" ,"Organism" ,"Length","Annotation","Cross.reference..DisProt.")
pbraUncharProteins <- pbraUniprot[pbraUniprot$Protein.names=="Uncharacterized protein",]
# almost all NA genes are uncharacterized. Good, can remove those.
pbra_allHighDomain_nm_rmUnchar <- pbra_allHighDomain_nm[!(pbra_allHighDomain_nm$target %in% pbraUncharProteins$yourlist),]
pbra_allHighDomain_rmUnchar <- pbra_allHighDomain[!(pbra_allHighDomain$Gene %in% pbraUncharProteins$yourlist),]
dim(pbra_allHighDomain_rmUnchar) # these are vep results without gconvert
dim(pbra_allHighDomain_nm_rmUnchar) # this has names from gconvert
length(unique(pbra_allHighDomain_rmUnchar$Gene)) # 1503  unique genes once unchar removed

## elut uniprot:
# do this once: (20171228)
write_clip(elut_allHighDomain_nm[elut_allHighDomain_nm$name=="N/A",]$alias)
# then go to uniprot:
# http://www.uniprot.org/uniprot/ search and save results to
# Old: 

elutUniprot <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/uniprotVEP_onlineSearch_20180207/elut.VEP_NAGenes.UniprotResults.50_250DP.20180207",header=T,sep="\t")
head(elutUniprot)
colnames(elutUniprot) <- c("yourlist","isomap","Entry", "Entry.name","Status","Protein.names" ,"Gene.names" ,"Organism" ,"Length","Annotation","Cross.reference..DisProt.")
elutUncharProteins <- elutUniprot[elutUniprot$Protein.names=="Uncharacterized protein",]
# almost all NA genes are uncharacterized. Good, can remove those.
elut_allHighDomain_nm_rmUnchar <- elut_allHighDomain_nm[!(elut_allHighDomain_nm$target %in% elutUncharProteins$yourlist),]
elut_allHighDomain_rmUnchar <- elut_allHighDomain[!(elut_allHighDomain$Gene %in% elutUncharProteins$yourlist),]
dim(elut_allHighDomain_rmUnchar) # 3238 these are vep results without gconvert
dim(elut_allHighDomain_nm_rmUnchar) # this has names from gconvert
length(unique(elut_allHighDomain_nm_rmUnchar$alias)) # 1619

############################## INTERSECTION ####################
# intersect the lists that have removed uncharacterized proteins already;
elut_pbra_vep_rmUnchar <- intersect(pbra_allHighDomain_rmUnchar$Gene,elut_allHighDomain_rmUnchar$Gene)
length(elut_pbra_vep_rmUnchar) # 900 genes intersect under 50/250
length(intersect(pbra_allHighDomain_nm$alias,elut_allHighDomain_nm$alias)) # 1648 intersect before unchar removed
elut_pbra_vep_nm_rmUnchar <- gconvert(elut_pbra_vep_rmUnchar,organism = "mfuro")

################ Combine elut and pbra to get overall gene info for later #######
bothspp_AllHighDomain <- rbind(elut_allHighDomain_rmUnchar,pbra_allHighDomain_rmUnchar)
####################### Gprofiler on VEP genes #####################
#gProfileR::set_base_url("https://biit.cs.ut.ee/gprofiler/") # this is the slightly older version (not as buggy)

########### THIS IS BEST SOLUTION: Remove uncharacterized proteins but KEEP IN those that have functions but not a gene symbol (e.g. Olfactory Receptors!) ######### 
#### Function to convert gProfiler results to GEM format for cytoscape
convertToGEM <- function(gProfilerResults){
  gProfilerResults$FDR <- gProfilerResults$p.value
  gProfilerResults$Phenotype <- 1
  gProfilerResults_forCytoscape <- gProfilerResults[,c("term.id","term.name","p.value","FDR","Phenotype","intersection")]
  colnames(gProfilerResults_forCytoscape) <- c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")
  return(gProfilerResults_forCytoscape)
}
#### ELUT + PBRA Intersect genes (removed uncharacterized proteins)
# This is for cytoscape:
g00a <- gprofiler(unique(as.character(elut_pbra_vep_nm_rmUnchar$target)),organism="mfuro",min_set_size = 3,max_set_size = 1000,min_isect_size = 2,domain= "known",correction_method = "gSCS",include_graph = T,src_filter = c("GO:BP","KEGG","REACTOME"),hier_filtering = "moderate")
# restrict to only BP, kegg reactome for cytoscape
# make it tab delimited gem format for cytoscape
# gem output is:
#GO.ID	Description	p.Val	FDR	Phenotype	Genes
#GO:0008150	biological_process	3.41e-118	3.41e-118	1	ENSMPUG00000000007,ENSMPUG0
# make dummy columns:
#g00a$FDR <- g00a$p.value
#g00a$Phenotype <- 1
#g00a_forCytoscape <- g00a[,c("term.id","term.name","p.value","FDR","Phenotype","intersection")]
#colnames(g00a_forCytoscape) <- c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")

g00a_forCytoscape <- convertToGEM(g00a)
## elut only: 
g00e <- gprofiler(unique(as.character(elut_allHighDomain_nm_rmUnchar$target)),organism="mfuro",min_set_size = 3,max_set_size = 1000,min_isect_size = 2,domain= "known",correction_method = "gSCS",include_graph = T,src_filter = c("GO:BP","KEGG","REACTOME"),hier_filtering = "moderate"); View(g00e)
g00e$FDR <- g00e$p.value
g00e$Phenotype <- 1
g00e_forCytoscape <- g00e[,c("term.id","term.name","p.value","FDR","Phenotype","intersection")]
colnames(g00e_forCytoscape) <- c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")




### pbra all vep genes: 
g00p <- gprofiler(unique(as.character(pbra_allHighDomain_nm_rmUnchar$target)),organism="mfuro",min_set_size = 3,max_set_size = 1000,min_isect_size = 2,domain= "known",correction_method = "gSCS",include_graph = T,src_filter = c("GO:BP","KEGG","REACTOME"),hier_filtering = "moderate"); View(g00p)
g00p$FDR <- g00p$p.value
g00p$Phenotype <- 1
g00p_forCytoscape <- g00p[,c("term.id","term.name","p.value","FDR","Phenotype","intersection")]
colnames(g00p_forCytoscape) <- c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")

# get way more hits if you use Gene Symbols. Losing N/A genes... 

datadir=paste("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/VEP_Results_CompareElutPbra/gProfiler_VEP_Results_50_250DPFilter/includesIndels")
dir.create(path=datadir,recursive = T)

write.table(g00e,paste(datadir,"/gProfiler.elutHighImpact.SnpsIndels.IncludesNA.ExcludesUncharProteinsFromUniprot.mfuro.50_250DPFilter.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(g00p,paste(datadir,"/gProfiler.pbraHighImpact.SnpsIndels.IncludesNA.ExcludesUncharProteinsFromUniprot.mfuro.50_250DPFilter.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(g00a,paste(datadir,"/gProfiler.elutAndpbraIntersect.HighImpact.SnpsIndels.IncludesNA.ExcludesUncharProteinsFromUniprot.mfuro.50_250DPFilter.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(g00a_forCytoscape,paste(datadir,"/gProfiler.elutAndpbraIntersect.HighImpact.SnpsIndels.IncludesNA.ExcludesUncharProteinsFromUniprot.mfuro.GEMFormat.ForCytoscape.HierModerate.50_250DPFilter.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(g00e_forCytoscape,paste(datadir,"/gProfiler.elut.AllInlcudingShared.VEP.HighImpact.SnpsIndels.IncludesNA.ExcludesUncharProteinsFromUniprot.mfuro.GEMFormat.ForCytoscape.HierModerate.50_250DPFilter.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(g00p_forCytoscape,paste(datadir,"/gProfiler.pbra.AllIncludingShared.VEP.HighImpact.SnpsIndels.IncludesNA.ExcludesUncharProteinsFromUniprot.mfuro.GEMFormat.ForCytoscape.HierModerate.50_250DPFilter.txt",sep=""),row.names=F,quote=F,sep="\t")
###################### Plot results of gProfiler ##################
# PLOT IN CYTOSCAPE (outside of R)
######## GO SUMMMARIES (THIS DOESN'T WORK YET: ###########
#require(GOsummaries)
#genes <- c(as.character((elut_pbra_vep_nm_rmUnchar[elut_pbra_vep_nm_rmUnchar$name!="N/A",]$target)))
#genes1 =  c("MAP2K3" ,  "SLC5A10",  "ERFE"   ,  "SMS"   ,   "SNED1"  ,  "OR5AU1"  , "CASP10"  , "SMAD3"  ,  "GLYATL2" , "IRAK4"   , "TMPRSS7" , "ABHD10" ,  "RPAP3",    "RAPH1"   , "HDAC7",    "SERPINB9", "TBX19" ,   "FARS2" ,   "ANKDD1A" , "RCSD1"   )
#genes2 = c("201890_at", "202503_s_at", "204170_s_at", "201291_s_at", 
           #"202589_at", "218499_at", "209773_s_at", "204026_s_at")
#gl <- list(List1=genes)
#gs_pbra <- gosummaries(x=gl,ordered_query = F,domain="known",organism="mfuro")
#plot(gs_pbra)
# is it just no results?
# ist here a 15 gene minimum?
######### GOSummaries DOESNT WORK
## okay want to add the gprofiler info to the vep info
# start with test and then do lapply or something 
#################### Unnest intersection of Gprofiler results: g01a (ENSMPUG) ########
library(tidyr)
library(dplyr)
######## elut and pbra intersection (g00a) #######
g00a2 <- g00a %>% 
  mutate(intersection = strsplit(as.character(intersection), ",")) %>% 
  unnest(intersection) # this makes one row for each gene in an intersection! Nice!
# Then can merge it with VEP info:
g00a3 <- merge(g00a2,bothspp_AllHighDomain,by.x="intersection",by.y="Gene")
# this is so cool! 
# want to condense consequences:
g00a3$plotConsequence <- "4. Other"
# if it has a stopgained (or more)
# note! grepl does indexing. you can do it either way with grep or grepl, but grepl is what works for frameshift AND stop gained 
g00a3[grepl("stop_gained",g00a3$Consequence),]$plotConsequence <- "2. Stop-Gained(+)"
g00a3[grepl("frameshift",g00a3$Consequence),]$plotConsequence <- "1. Frameshift(+)"
g00a3[grepl("frameshift",g00a3$Consequence) & grepl("stop_gained",g00a3$Consequence),]$plotConsequence <- "3. Frameshift+Stop-Gained"
# others:
#[1] "start_lost"                                
#[2] "splice_donor_variant,missense_variant"     
#[3] "splice_donor_variant,inframe_insertion"    
#[4] "splice_acceptor_variant,synonymous_variant"
#[5] "stop_lost,inframe_deletion"
library(stringr)
# need to wrap labels:
g00a3$term.name.wrapped = str_wrap(g00a3$term.name, width = 50)

## also want to collapse some terms if they share >95% of same genes:

###################### FOR MANUSCRIPT: Compare VEP and Pseudopipe ##########
#### read in pseudopipe results: (not changed by 50/150 or 50/250 filters bc from de novo genome)####
# 20170102: updated it to remove frac > 0.75 filter (changed it in the USETHIS_reorderSteps_newSetDiff script; new files have "noFracFilter" in name)
elut_only_finalPgenes_gconvert <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/pseudoPipe_20170307/pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_only/elut_only_finalPgenes_gconvert_NonZeroDisruptions.noFracFilter.txt",sep="\t",quote="",header=T)
elut_only_finalPgenes_gconvert$spp <- "elut"
pbra_only_finalPgenes_gconvert <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/pseudoPipe_20170307/pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/pbra_only/pbra_only_finalPgenes_gconvert_NonZeroDisruptions.noFracFilter.txt",quote="",header=T,sep="\t")
pbra_only_finalPgenes_gconvert$spp <- "pbra"
# shared:
ep_Pgenes_ep_FINAL_gconvert<- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/pseudoPipe_20170307/pseudoPipeResults/ResultTablesAfterFiltering_ElutPbra/elut_and_pbra/elut_and_pbra_shared_finalPgenes_gconvert_NonZeroDisruptions.noFracFilter.txt",quote="",header=T,sep="\t")
ep_Pgenes_ep_FINAL_gconvert$spp <- "shared"
elutNames <- names(ep_Pgenes_ep_FINAL_gconvert)[1:18]
pbraNames <- c(names(ep_Pgenes_ep_FINAL_gconvert)[1:4],names(ep_Pgenes_ep_FINAL_gconvert)[19:32])
colnames(elut_only_finalPgenes_gconvert) <- c(elutNames,"spp")
colnames(pbra_only_finalPgenes_gconvert) <- c(pbraNames,"spp")

allElutPgenes <- rbind(elut_only_finalPgenes_gconvert,ep_Pgenes_ep_FINAL_gconvert[,c(elutNames,"spp")])
allPbraPgenes <- rbind(pbra_only_finalPgenes_gconvert,ep_Pgenes_ep_FINAL_gconvert[,c(pbraNames,"spp")])

# number of elut only unique genes: (includes NA genes )
length(unique(elut_only_finalPgenes_gconvert$alias)) # 467 w/out frac filter (563 with frac 0.75 filter -- increases because more are moved to elut-only)
length(unique(pbra_only_finalPgenes_gconvert$alias)) # 538 /wout frac filter; (328 with frac filter)
length(unique(ep_Pgenes_ep_FINAL_gconvert$alias)) # 725 (232 w/ frac 0.75 filter)
############# Lookup Pseudopipe N/As in Uniprot and remove unchar proteins ########
########## ELUT
#write_clip(allElutPgenes[allElutPgenes$name=="N/A",]$target) # 249 N/A elut genes (115 with frac  0.75 filter) # DO ONCE
# go to Uniprot, get identifiers, download results:
# read in results: 
allElutPgenes_NA_UniprotResults <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/Uniprot_Pseudopipe_NA_Results/withoutFracFilter/elut.noFracFilter.Pseudopipe_all_NAGenes.UniprotResults.20180102",header=T,sep="\t")
colnames(allElutPgenes_NA_UniprotResults) <- c("yourlist","isomap","Entry", "Entry.name","Status","Protein.names" ,"Gene.names" ,"Organism" ,"Length","Annotation","Cross.reference..DisProt.")
# pull out uncharacterized proteins
allElutPgenes_NA_UniprotResults_UncharProteins <- allElutPgenes_NA_UniprotResults[allElutPgenes_NA_UniprotResults$Protein.names=="Uncharacterized protein",]
dim(allElutPgenes_NA_UniprotResults_UncharProteins) # 133/249 are unchar
# use this downstream:
allElutPgenes_rmUnchar <- allElutPgenes[!(allElutPgenes$target %in% allElutPgenes_NA_UniprotResults_UncharProteins$yourlist),]
length(unique(allElutPgenes$target)) # 1192 
length(unique(allElutPgenes_rmUnchar$target)) # 1060


#### PBRA
#length(allPbraPgenes[allPbraPgenes$name=="N/A",]$target) # 260
#write_clip(allPbraPgenes[allPbraPgenes$name=="N/A",]$target) 
# go to Uniprot, get identifiers, download results:
# read in results:
allPbraPgenes_NA_UniprotResults <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/Uniprot_Pseudopipe_NA_Results/withoutFracFilter/pbra.noFracFilter.Pseudopipe_all_NAGenes.UniprotResults.20180102",header=T,sep="\t")
colnames(allPbraPgenes_NA_UniprotResults) <- c("yourlist","isomap","Entry", "Entry.name","Status","Protein.names" ,"Gene.names" ,"Organism" ,"Length","Annotation","Cross.reference..DisProt.")
# pull out uncharacterized proteins
allPbraPgenes_NA_UniprotResults_UncharProteins <- allPbraPgenes_NA_UniprotResults[allPbraPgenes_NA_UniprotResults$Protein.names=="Uncharacterized protein",]
dim(allPbraPgenes_NA_UniprotResults_UncharProteins) # 132/260 are unchar
# use this downstream:
allPbraPgenes_rmUnchar <- allPbraPgenes[!(allPbraPgenes$target %in% allPbraPgenes_NA_UniprotResults_UncharProteins$yourlist),]
dim(allPbraPgenes) # 2550 (up from 590 when frac filter applied)
dim(allPbraPgenes_rmUnchar) # 2318  (up from 510)
######### final counts of unique pseudopipe genes once UncharProt removed ######
length(unique(allPbraPgenes_rmUnchar$target))
length(unique(allElutPgenes_rmUnchar$target))
length(intersect(allElutPgenes_rmUnchar$target,allPbraPgenes_rmUnchar$target))
######## Compare VEP and Pseudopipe: #####
elut_allHighDomain_rmUnchar
pbra_allHighDomain_rmUnchar
#### USE THESE downstream
# need to put all elut Pgenes together
############# Merge VEP and Pseudopipe (post removing unchar proteins) #########
elutVEP_Ppipe_Merge <- (merge(allElutPgenes_rmUnchar,elut_allHighDomain_rmUnchar,by.x="target",by.y="Gene"))
pbraVEP_Ppipe_Merge <-  (merge(allPbraPgenes_rmUnchar,pbra_allHighDomain_rmUnchar,by.x="target",by.y="Gene"))
sharedElutPbraVEP_Pipe_Merge <- intersect(elutVEP_Ppipe_Merge$target,pbraVEP_Ppipe_Merge$target)
length(unique(elutVEP_Ppipe_Merge$target)) # new: 113
length(unique(pbraVEP_Ppipe_Merge$target)) # same before and after filter 133
length(unique(sharedElutPbraVEP_Pipe_Merge)) # 45

############### gProfiler 01: on VEP/Pseudopipe Intersection ###########
datadir2=paste("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/VEP_Results_CompareElutPbra/gProfiler_VEP_Results_50_250DPFilter/includesIndels/")
dir.create(path=datadir2,recursive = T)

### these are shared between elut and pbra AND found in VEP and Ppipe: no results! 
g01a <-  gprofiler(unique(as.character(sharedElutPbraVEP_Pipe_Merge)),organism="mfuro",min_set_size = 3,max_set_size = 1000,min_isect_size = 2,domain= "known",correction_method = "gSCS",include_graph = T,src_filter = c("GO:BP","KEGG","REACTOME")); View(g01a)

##### elut (all VEP + pseudopipe pbra genes; not just pbra-only)

g01e <-  gprofiler(unique(as.character(elutVEP_Ppipe_Merge$target)),organism="mfuro",min_set_size = 3,max_set_size = 1000,min_isect_size = 2,domain= "known",correction_method = "gSCS",include_graph = T,src_filter = c("GO:BP","KEGG","REACTOME"),hier_filtering = "moderate"); View(g01e) # not doing hierarchical filtering
g01e$FDR <- g01e$p.value
g01e$Phenotype <- 1
g01e$spp <- "sea otter"
g01e_forCytoscape <- g01e[,c("term.id","term.name","p.value","FDR","Phenotype","intersection")]
colnames(g01e_forCytoscape) <- c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")
##### pbra (all VEP + pseudopipe pbra genes; not just pbra-only)
g01p <-  gprofiler(unique(as.character(pbraVEP_Ppipe_Merge$target)),organism="mfuro",min_set_size = 3,max_set_size = 1000,min_isect_size = 2,domain= "known",correction_method = "gSCS",include_graph = T,src_filter = c("GO:BP","KEGG","REACTOME"),hier_filtering = "moderate"); View(g01p)
g01p$FDR <- g01p$p.value
g01p$Phenotype <- 1
g01p$spp <- "giant otter"
g01p_forCytoscape <- g01p[,c("term.id","term.name","p.value","FDR","Phenotype","intersection")]
colnames(g01p_forCytoscape) <- c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")

write.table(elutVEP_Ppipe_Merge,paste(datadir2,"Elut.Vep.Pseudopipe.NoFracFilter.Intersect.USETHIS.AllInfo.ForSupplement.50_250DPFilter.txt",sep=""), row.names = F,quote=F,sep="\t")
write.table(pbraVEP_Ppipe_Merge,paste(datadir2,"Pbra.Vep.Pseudopipe.NoFracFilter.Intersect.USETHIS.AllInfo.ForSupplement.50_250DPFilter.txt",sep=""), row.names = F,quote=F,sep="\t")

write.table(g01p_forCytoscape,paste(datadir2,"Pbra.Vep.NoFracFilter.Pseudopipe.SharedGenes.gProfiler.HierModerate.forCytoscape.GEM.50_250DPFilter.txt",sep=""), row.names = F,quote=F,sep="\t")
write.table(g01e_forCytoscape,paste(datadir2,"Elut.Vep.NoFracFilter.Pseudopipe.SharedGenes.gProfiler.HierModerate.forCytoscape.GEM.50_250DPFilter.txt",sep=""), row.names = F,quote=F,sep="\t")

######### Look at sensory genes: ##########
# get gene names: 
require(dplyr)
require(tidyr)
g01e2 <- g01e %>% 
  mutate(intersection = strsplit(as.character(intersection), ",")) %>% 
  unnest(intersection) # this makes one row for each gene in an intersection! Nice!

g01p2 <- g01p %>% 
  mutate(intersection = strsplit(as.character(intersection), ",")) %>% 
  unnest(intersection) # this makes one row for each gene in an intersection! Nice!
#### IMPORTANT FOR RESULTS: THESE ARE THE ELUT VEP/PPipe genes that make up the enrichments (p < 0.05): ########
ElutSensoryGenes <- elutVEP_Ppipe_Merge[elutVEP_Ppipe_Merge$target %in% g01e2$intersection,]
ElutSensoryGenes$name
PbraSensoryGenes <- pbraVEP_Ppipe_Merge[pbraVEP_Ppipe_Merge$target %in% g01p2$intersection,]
PbraSensoryGenes$name
## count distinct ENSMPUG genes:
length(unique(PbraSensoryGenes$target)) # 15
length(unique(ElutSensoryGenes$target)) # 15
#### SO : results didn't change after 50/250 filter! Still are 12 ORs + 3 other genes 
# look up the NAs in uniprot to check:
elutUniprot[elutUniprot$yourlist %in% ElutSensoryGenes$target,] # all are olfactory receptors 
pbraUniprot[pbraUniprot$yourlist %in% PbraSensoryGenes$target,] # all are olfactory receptors 
########## lowest pvalues: ##########
g01e[g01e$p.value==min(g01e$p.value),] # sea otter lowest p value
g01p[g01p$p.value==min(g01p$p.value),] # detection of chemical stimulus involved in sensory perception   
################# Do research on these sensory genes
