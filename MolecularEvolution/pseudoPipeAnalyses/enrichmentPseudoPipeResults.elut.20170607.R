### are pseudogenes enriched for anything?
 ### this is elut:elut
pseudo <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/pseudoPipe_20170307/pseudoPipeResults/pseudoPipe_elut_elutPep/enhydra_lutris_ppipe_elutPep.allPgenes.20170321.txt",header=T,quote="",comment="")
dim(pseudo)
head(pseudo)
unique(pseudo$type)
# FRAG DUP  PSSD -- most interested in PSSD
sum(pseudo$type=="PSSD")
sum(pseudo$type=="DUP")
sum(pseudo$type=="FRAG")
#### info on otter genes ####
# Info on all elut genes
allElutGenes_fullInfo <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Annotation/MakerRound5_blast1e-06_USETHIS_Nov11/FULLINFO_GradedGeneSets_Nov17/mRNA_FullInfo_withGrades.Nov17.txt",sep=";",header=T,quote="",comment.char="")
# do some extra formatting to it:
allElutGenes_fullInfo$jbrowse <- paste(allElutGenes_fullInfo$scaffold,":",allElutGenes_fullInfo$start,"-",allElutGenes_fullInfo$stop,sep="")
library(stringr)
head(allElutGenes_fullInfo)
int1e <- str_split_fixed(allElutGenes_fullInfo$blastInfo,":",2)
head(int1e)
allElutGenes_fullInfo$similarToSymbol <- str_split_fixed(int1e[,1],"to ",2)[,2]
allElutGenes_fullInfo$similarToSymbol
head(allElutGenes_fullInfo)

#### merge pseudo info:
pseudo_info <- merge(pseudo,allElutGenes_fullInfo,by.x="query",by.y="name")
head(pseudo_info)
View(pseudo_info[pseudo_info$type=="PSSD",])
# unique gene symbols associated with PSSD pseudogenes:
uniqPSSD <- unique(pseudo_info[pseudo_info$type=="PSSD",]$similarToSymbol)
# some of these gene symobls are kind of messed up -- fix/
### next try gprofiler or something with these? 
require(gProfileR)

  # genesMod1=list of significant genes
  # organism= ref genome
  # correction method="gSCS" or "fdr". Reimand et al 2007 suggest gSCS is better
  # domain size = "annotated" or "known"
  # custom_bg=background list (all genes in genome)
  # for all parameters
  # ?gProfileR
  
GO <- gprofiler(uniqPSSD, organism = "mustelaputoriusfuro", ordered_query = F, significant = T, exclude_iea = F, underrep = F, evcodes = F, region_query = F, max_p_value = 1, min_set_size = 3, min_isect_size = 3, correction_method = "gSCS", hier_filtering = "none", domain_size = "known", custom_bg = c(allElutGenes_fullInfo$similarToSymbol), numeric_ns = "", png_fn = T, include_graph = T, src_filter = NULL); View(GO)
