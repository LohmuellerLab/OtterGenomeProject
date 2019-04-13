######################## Initial Tests of Daub et al.'s method for gene set enrichment analysis ###############
## This is the latest script (20171219), where artifacts have been filtered out, and SO HAVE ~ 1400 genes flagged as having high impact VEP snps or indels based on alignment to ferret (in either elut or pbra)
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/CompGenomicsScripts/newPostCodemlAnalysis_DaubExcoffier_20170525/polysel/")
#### NOTE: the setObj and setInfo files will work for all three branches! Just need to redo the codeml stuff. 
branch="otters"
pamlDate="20170922"
rundate="20171018" # polysel run date
# paper: https://academic.oup.com/mbe/article/34/6/1391/3053366
# github: https://github.com/CMPG/polysel
#projectname=paste(branch,"paml",pamlDate,"polyselRun",rundate,sep=".")
TNamespp1="hsap"
TNamespp2="cfam" # species used for TName labels ("cfam","hsap","fcat" are options)
# things you need (generic for all runs:)
# 1. gene2ensembl: List of gene UIDs from NCBI and their corresponding ensembl IDs (what org to use? human?)
gene2ensembl <- read.table("inputFormatting/20170525/gene2ensembl",comment.char = "",header=T)
dim(gene2ensembl) # 650676      7
head(gene2ensembl)

# 2. biosystems_gene_all: Set of all biosystems and genes that are in them. Note that conserved biosystems include gene IDs from all species they know about
# only need to do this once
biosystems_gene_all <- read.table("inputFormatting/20170525/biosystems_gene_all",header=F)
colnames(biosystems_gene_all) <- c("bsid","geneUID","score")

# 3. conserved Biosystems: BSIDs for the conserved biosystems (or any biosystems you want, could be organism specific). (https://www.ncbi.nlm.nih.gov/biosystems -- search for type: pathway, and restrict to conserved biosystems) # a lot of these are from bacteria/plants/etc. but we will remove them through filtering. 
# conserved biosystem IDs only:
conserved_bsid <- read.table("inputFormatting/20170525/conservedBiosystems.fromWeb.txt")
head(conserved_bsid)
colnames(conserved_bsid) <- "bsid"
dim(conserved_bsid) # 32990
# conserved biosystem with names:
conserved_bsid_names <- read.csv("inputFormatting/20170525/conservedBiosystems.Info.20170525.csv",header=T)
dim(conserved_bsid_names) # 32990 
# sanity check: are these the same bsids?
sum(conserved_bsid_names$BSID %in% conserved_bsid$bsid) # 32990, good.

head(conserved_bsid)
colnames(conserved_bsid) <- "bsid"
dim(conserved_bsid) # 32990
# conserved biosystem with names:
conserved_bsid_names <- read.csv("inputFormatting/20170525/conservedBiosystems.Info.20170525.csv",header=T)
dim(conserved_bsid_names) # 32990 
# sanity check: are these the same bsids?
#sum(conserved_bsid_names$BSID %in% conserved_bsid$bsid) # 32990, good.

# use a priori traits:
### 20171219: when searching ncbi biosystems, am restricting search term to "Description" field; this helps reduce numbers and greatly improves specificity to the terms!!
# "Description" is way better
# Still need to restrict to "conserved", kegg/go, and remove plant, yeast, bacteria, etc.

osmo <- read.csv("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PolyselPathwayLitReview/PathwayUIDs/searchTermDescription/osmoregulation.Description.csv")
thermo <- read.csv("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PolyselPathwayLitReview/PathwayUIDs/searchTermDescription/thermoregulation.Description.csv")
hair.follicle <- read.csv("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PolyselPathwayLitReview/PathwayUIDs/searchTermDescription/hair.follicle.description.csv")
diet <- read.csv("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PolyselPathwayLitReview/PathwayUIDs/searchTermDescription/diet.Description.csv")
sensory <- read.csv("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PolyselPathwayLitReview/PathwayUIDs/searchTermDescription/sensory.perception.Description.csv")
hypoxia <- read.csv("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PolyselPathwayLitReview/PathwayUIDs/searchTermDescription/hypoxia.Description.csv")

traits <- rbind(thermo,osmo,hair.follicle,diet,sensory,hypoxia)
dim(traits) # 914

traits <- traits[traits$Type=="conserved biosystem",]
dim(traits)  # 178 conserved pathways (look at them)

# want to restrict to GO and KEGG only (BIOCYC is all plants/bacteria, not right)
traits <- traits[traits$Source=="GO" | traits$Source=="KEGG",]
dim(traits) # 105
# get rid of anything with plant or bacteria or yeast:
# plants:
traits <- traits[grep("plant",traits$Name,ignore.case=T,invert = T),]
dim(traits)
traits <- traits[grep("bact",traits$Name,ignore.case=T,invert = T),]
dim(traits)
traits <- traits[grep("coli",traits$Name,ignore.case=T,invert = T),]
dim(traits)
traits <- traits[grep("yeast",traits$Name,ignore.case=T,invert = T),]
dim(traits) # 105 <-- this is a good set! actually specific, not about plants, etc. 
#View(traits)


#################################### 20170705 otters paml results #####################################
# specific to this run: 
# 4. The results from PAML. These are branch specific for now (have to figure out fdr later)
# 20170922: I am starting with the otters branch from the gblocks filtering, with artifacts removed (inspected all genes with p < 0.01)
flag="otters"
#rep="1"
pamlDate="20170922"

# this has already excluded BEB > 1, length < 120bp, and artifacts that failed visual inspection
### THIS IS THE LATEST (20171228): N/As are included if they have a uniprot known function (olfactory receptors, transporters, things like that)
# artifacts (NA and regular) have been removed
# so have "VEP high impact" genes.
codeml <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/otters.fixedNAs.RemovedVEPHighImpactGenes.50_250DPFilter.putativelyReal.andUnexaminedGenes.20170922.forPolysel.USETHIS.txt",header=T,sep="\t") 

dim(codeml) # 11083

############################### Remove LGR4 Gene #######################################
dim(codeml)
codeml <- codeml[codeml$Symbol!="LGR4",]
dim(codeml) # make sure this went down by 1 (11083 --> 11082)
############################## GET NCBI IDs #########################
cfamtt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T)
hsaptt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T,quote="")


################################# 0. PREPARE INPUT #############################
#### 0a. want to get UIDs for all genes: Human_Gene ID has a few mistakes/ambiguities
# instead going to try human protein ID, which is also in the gene2ensembl thing. 

# get rid of unnamed proteins:
gene2ensembl_nounnamed <- gene2ensembl[gene2ensembl$Ensembl_protein_identifier!="-",]

dim(gene2ensembl)
dim(gene2ensembl_nounnamed) # first merge based on the protein IDs and gene IDs together
# then there are some that don't have protein ID but do have gene ID
head(codeml)
codeml$TName.No.1 <- gsub(perl=T,"\\.[0-9]+","",codeml$TName)
head(codeml) # got rid of .1 digit

### this should work for any codeml run:
# a. first get the ones that match based on TName:
# get gene IDs for cfam and hsap:
codeml_hsap <- merge(codeml, unique(hsaptt),by.x = "TName",by.y="Transcript",all=F)
dim(codeml_hsap) # 12570 >> 10969 post VEP 
codeml_cfam <- merge(codeml, unique(cfamtt),by.x="TName",by.y="cfam",all=F)
dim(codeml_cfam) # 359 >> 97 post VEP 
dim(codeml_cfam) + dim(codeml_hsap) # 11066 (total, matches original codeml dim, good)

# b. get geneIDs based on gene name (note, same UID for multiple isoforms)
codeml_hsap2 <- merge(codeml_hsap, unique(gene2ensembl_nounnamed[,c("GeneID","Ensembl_gene_identifier")]),by.x="GeneNo.1",by.y="Ensembl_gene_identifier",all=F)
dim(codeml_hsap2) # 12539 # lost 24 genes -- why?
# merge leftovers using protein id: 
codeml_hsap3 <- merge(codeml_hsap[!codeml_hsap$TName %in% codeml_hsap2$TName,], unique(gene2ensembl_nounnamed[,c("GeneID","Ensembl_protein_identifier")]),by.x="Protein",by.y="Ensembl_protein_identifier",all=F)
dim(codeml_hsap3) # gain 2 more genes
## combine:
codeml_hsap_combo <- rbind(codeml_hsap2,codeml_hsap3)

dim(codeml_hsap3) # gain 2 more genes
codeml_cfam2 <-  merge(codeml_cfam, unique(gene2ensembl_nounnamed[,c("GeneID","Ensembl_gene_identifier")]),by.x="Dog_GeneNo.1",by.y="Ensembl_gene_identifier",all=F)
dim(codeml_cfam2) # 278 # lost 81 genes
codeml_cfam3 <- merge(codeml_cfam[!codeml_cfam$TName %in% codeml_cfam2$TName,], unique(gene2ensembl_nounnamed[,c("GeneID","Ensembl_protein_identifier")]),by.x="Dog_Protein",by.y="Ensembl_protein_identifier",all=F)
dim(codeml_cfam3) # gain 0 more genes

codeml_cfam_combo <- rbind(codeml_cfam2,codeml_cfam3)
# caculate genes lost:
dim(codeml)[1] - (dim(codeml_cfam_combo)[1] + dim(codeml_hsap_combo)[1]) # 32 genes lost (not bad) (previously was 108 genes lost)

# combine hsap and cfam:
codeml_all <- rbind(codeml_cfam_combo[,!names(codeml_cfam_combo) %in% c("Dog_Gene","Dog_ProteinNo.1","Symbol.y","Dog_Protein","TName.No.1","Dog_GeneNo.1")],codeml_hsap_combo[,!names(codeml_hsap_combo) %in% c("GeneNo.1","Symbol.y","Protein","TName.No.1","ProteinNo.1","Gene")])
dim(codeml_all)           

#### 0b. Prepare paml results (get lnL4) NOW USING maxLRT (3 reps) #### 

# get delta lnL4 (fourth root of delta lnL)
codeml_all$dlnl4 <- (codeml_all$maxLRT)^(1/4)
head(codeml)


######## add in extra info to codeml_uid:

#### Od. Make OBJINFO file with  ‘objID’, ‘objStat’ and ‘objName’ (and then a couple other optional files, I want group, cluster,TName, bebCount, etc. so that I can keep track of everything too) #### 
objInfo <- codeml_all[,c("GeneID","dlnl4","Symbol.x")] # Symbol.x may end up being Symbol
head(objInfo)
colnames(objInfo) <- c("objID","objStat","objName")
head(objInfo)

############ 0d.2 add in bonus columns of extra gene info ################
objInfo <- cbind(objInfo,codeml_all[,c("group","cluster","branch.x","replicate","sppCount",paste(flag,"treeLength",sep="."),paste(flag,"bebCount",sep="."),"TName","seqLength_mayInclGaps","significant")])

# for now not adding : ,"scaffold","start","stop","strand")])
head(objInfo)

## make a folder for this project in the polysel/data project 
datadir=paste("data/otters/paml.",pamlDate,"-traits-postVEP_50_250DPFilter-WithoutLGR4/",branch,sep="")
dir.create(path=datadir,recursive = T)

write.table(objInfo,paste(datadir,"/ObjInfo.txt",sep=""),row.names=F,quote=F,sep="\t")

#### 0e. Make SETINFO file with SETID and Setname (download the csv file from ncbi biosystems) ####
setInfo <- traits[,c("BSID","Name","Source")]

colnames(setInfo) <- c("setID","setName","setSource")
write.table(setInfo,paste(datadir,"/SetInfo.txt",sep=""),row.names=F,quote=F,sep="\t")


#### Of. make SETOBJ file with BSID and gene ID (only conserved biosystems) #### 
head(conserved_bsid)
biosystems_gene_conserved <- biosystems_gene_all[(biosystems_gene_all$bsid %in% traits$BSID),]

dim(biosystems_gene_conserved)
head(biosystems_gene_conserved)
dim(biosystems_gene_conserved) # 18821542
dim(conserved_bsid) # 32990
length(unique(biosystems_gene_conserved$bsid)) # 23395 so 9595 are missing
head(biosystems_gene_conserved)

## make SETOBJ file:
setObj <- biosystems_gene_conserved[,c("bsid","geneUID")]
colnames(setObj) <- c("setID","objID")
# this includes many pathways that don't have any of my genes in them, that's okay will be taken care of by the scripts. 
write.table(setObj,paste(datadir,"/SetObj.txt",sep=""),row.names=F,quote=F,sep="\t")

