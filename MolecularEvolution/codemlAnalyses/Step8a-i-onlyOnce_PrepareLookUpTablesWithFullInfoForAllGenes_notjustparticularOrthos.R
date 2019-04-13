########## PREPARE LOOK UP TABLES FOR ALL ORTHOLOG CLUSTERS ################
cfamt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.lookuptable.txt",stringsAsFactors = F)
colnames(cfamt) <- c("Protein","Gene","Transcript")
cfamt$ProteinNo.1 <- gsub("\\.[0-9]+","",cfamt$Protein)
cfamt$GeneNo.1 <- gsub("\\.[0-9]+","",cfamt$Gene)
colnames(cfamt) <- c("Dog_Protein","Dog_Gene","cfam","Dog_ProteinNo.1","Dog_GeneNo.1")


hsapt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.lookuptable.txt",stringsAsFactors = F)
colnames(hsapt) <- c("Protein","Gene","Transcript")
hsapt$ProteinNo.1 <- gsub("\\.[0-9]+","",hsapt$Protein)
hsapt$GeneNo.1 <- gsub("\\.[0-9]+","",hsapt$Gene)

write.table(hsapt$GeneNo.1,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.GeneNames.FORGCONVERT.20170702.txt",quote=F,row.names=F,col.names=F)
write.table(cfamt$Dog_GeneNo.1,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.GeneNames.FORGCONVERT.20170702.txt",quote=F,row.names=F,col.names=F)

#### THEN USE THESE IN gCONVERT TO GET SYMBOLS FOR ALL GENE NAMES (then link back to transcripts)
## 

## then read them back in and link up with transcripts and protein names to have full info (can then merge with codeml results)

hsap_wSymbols <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.AllGenes.WithSymbols.FromGConvert.20170702.txt",header = F)
colnames(hsap_wSymbols) <- c("GeneNo.1","Symbol")
head(hsap_wSymbols)
dim(hsap_wSymbols) #102915
hsapt_wSymbols <- merge(hsapt,hsap_wSymbols,by.x="GeneNo.1",by.y="GeneNo.1")
head(hsapt_wSymbols)


cfam_wSymbols <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.AllGenes.WithSymbols.FromGConvert.20170702.txt",header = F)
colnames(cfam_wSymbols) <- c("GeneNo.1","Symbol")
head(cfam_wSymbols)
dim(cfam_wSymbols) #25157
cfamt_wSymbols <- merge(cfamt,cfam_wSymbols,by.x="Dog_GeneNo.1",by.y="GeneNo.1")
head(cfamt_wSymbols)

write.table(cfamt_wSymbols,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",quote=F,row.names=F,col.names=T)
write.table(hsapt_wSymbols,"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",quote=F,row.names=F,col.names=T)
