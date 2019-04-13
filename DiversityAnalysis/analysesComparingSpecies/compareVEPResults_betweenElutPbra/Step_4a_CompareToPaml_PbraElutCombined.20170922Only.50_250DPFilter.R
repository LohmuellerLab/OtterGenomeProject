### Script
# Want to find VEP results and see if any of those genes intersect with my pos selection results
# if so, mark as "VEP" and remove. how to deal with? need symbol 
# 
######### PORTHO RESULTS AS LOOKUP TABLE ##########
portho <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/proteinOrtho_20170517_newspp/15spp_POrtho_LongestIsosSelected_pruneTree_20170601/myproject.fixedHeader.poff",header=T,na.strings = c("*",","),stringsAsFactors = F)

# also need hsaptt and cfamtt lookup tables:
# ths includes seqLength, treelength, etc. 
cfamtt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/cfam.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T)
hsaptt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/hsap.FINAL.lookuptable.InclSymbols.20170702.USETHIS.txt",header=T,quote="")
mfurtt <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/Homologs_CDS_GFF_ALLSPP_Oct26_2016/CDS_GFF_fromEnsembl_forPORTHO_Oct26/lookupTables/mfur.lookuptable.txt",header=F,quote="")
colnames(mfurtt) <- c("Protein","Gene","Transcript")
# get rid of 1s: 
mfurtt$Protein.No.1 <- gsub(".[0-9]$","",mfurtt$Protein,perl=T)
mfurtt$Transcript.No.1 <- gsub(".[0-9]$","",mfurtt$Transcript,perl=T)
mfurtt$Gene.No.1 <- gsub(".[0-9]$","",mfurtt$Gene,perl=T)


########## VEP RESULTS  ##########
# high impact, must have an associated protein domain: 
### ELUT: 

vep_indels_high_domain_elut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/01_Elut_CA_Gidget.raw_variants.cdsSequence.AllScaffsConcat.20171006.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
colnames(vep_indels_high_domain_elut) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(vep_indels_high_domain_elut) # 5167 indels


vep_snps_high_domain_elut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/01_Elut_CA_Gidget.raw_variants.cdsSequence.AllScaffsConcat.20171006.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
colnames(vep_snps_high_domain_elut) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(vep_snps_high_domain_elut) # 711 (was 679)

allHighDomain_elut <- rbind(vep_indels_high_domain_elut,vep_snps_high_domain_elut)

#### PBRA


vep_indels_high_domain_pbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/PteBra_1.raw_variants.cdsSequence.AllScaffsConcat.20171206.HQsites.Only.i-4.PASS-ONLY.HF.Indels.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
colnames(vep_indels_high_domain_pbra) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(vep_indels_high_domain_pbra) # 4544


vep_snps_high_domain_pbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/filteredVepOutput1_HighImpact/PteBra_1.raw_variants.cdsSequence.AllScaffsConcat.20171206.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.HomozygousAlt.VEP.output.minimalInfo.HIGH.DOMAIN.CANON.tbl",stringsAsFactors = F)
colnames(vep_snps_high_domain_pbra) <- c("Uploaded_variation"	,"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence"	,"cDNA_position"	,"CDS_position"	,"Protein_position",	"Amino_acids"	,"Codons",	"Existing_variation",	"Extra")
dim(vep_snps_high_domain_pbra) # 705

allHighDomain_pbra <- rbind(vep_indels_high_domain_pbra,vep_snps_high_domain_pbra)
dim(allHighDomain_pbra) 
# look at overlap with positive selection results
########### PAML RESULTS (20170922 for polysel; doing 20170717 in a separate script)  #########
codemlElut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_elut/preVEP/elut.putativelyReal.andUnexaminedGenes.20170922.forPolysel.txt",header=T,quote="",sep="\t",stringsAsFactors = F) # gblocks

codemlOtters <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/preVEP/otters.putativelyReal.andUnexaminedGenes.20170922.forPolysel.txt",header=T,quote="",sep="\t",stringsAsFactors = F)


############# have to think about what to do here. Where is original file? 
codemlPbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/preVEP/pbra.putativelyReal.andUnexaminedGenes.20170922.forPolysel.txt",header=T,quote="",sep="\t",stringsAsFactors = F)

################# Get ferret ID for each paml cluster ##########
######## THE ISSUE
# Vep gives you results in terms of ferret ID
# Paml I used human and dog as IDs
# need to go Human/Dog TName (codeml) >> Human/Dog PName (tt) >> Mfur Pname (portho) >> Mfur GName (mfurtt)
# merge things: 
# portho and mfur look up table
mer1 <- merge(portho,mfurtt,by.x="mfur.faa",by.y="Protein.No.1",all.x=F,all.y=F) # this combines portho and ferret info 
# codeml and hsap/cfam look up tables
mer2_hsap_elut <- merge(codemlElut,unique(hsaptt),by.x="TName",by.y="Transcript",all.x = F,all.y=F)
mer2_cfam_elut <- merge(codemlElut,unique(cfamtt),by.x="TName",by.y="cfam",all.x=F,all.y=F)
mer2_hsap_otters <- merge(codemlOtters,unique(hsaptt),by.x="TName",by.y="Transcript",all.x = F,all.y=F)
mer2_cfam_otters <- merge(codemlOtters,unique(cfamtt),by.x="TName",by.y="cfam",all.x=F,all.y=F)
mer2_hsap_pbra <- merge(codemlPbra,unique(hsaptt),by.x="TName",by.y="Transcript",all.x = F,all.y=F)
mer2_cfam_pbra <- merge(codemlPbra,unique(cfamtt),by.x="TName",by.y="cfam",all.x=F,all.y=F)
# these should add up to dim(codemlElut_20170922) YES THEY DO! :)

## adds in ferret protein info to codeml results
mer3_hsap_elut <- merge(mer2_hsap_elut,portho[,c("hsap.faa","mfur.faa")],by.x="ProteinNo.1",by.y="hsap.faa") ## adds in ferret protein info to codeml results
mer3_cfam_elut <- merge(mer2_cfam_elut,portho[,c("cfam.faa","mfur.faa")],by.x="Dog_ProteinNo.1",by.y="cfam.faa") 

mer3_hsap_otters <- merge(mer2_hsap_otters,portho[,c("hsap.faa","mfur.faa")],by.x="ProteinNo.1",by.y="hsap.faa") ## adds in ferret protein info to codeml results
mer3_cfam_otters <- merge(mer2_cfam_otters,portho[,c("cfam.faa","mfur.faa")],by.x="Dog_ProteinNo.1",by.y="cfam.faa") 

mer3_hsap_pbra <- merge(mer2_hsap_pbra,portho[,c("hsap.faa","mfur.faa")],by.x="ProteinNo.1",by.y="hsap.faa") ## adds in ferret protein info to codeml results
mer3_cfam_pbra <- merge(mer2_cfam_pbra,portho[,c("cfam.faa","mfur.faa")],by.x="Dog_ProteinNo.1",by.y="cfam.faa") 

## adds ferret transcript and gene into to codeml results:
mer4_hsap_elut <- merge(mer3_hsap_elut,mfurtt[,c("Protein.No.1","Transcript.No.1","Gene.No.1")],by.x="mfur.faa",by.y="Protein.No.1")
dim(mer4_hsap_elut) # not all of them have ferret orthologs
mer4_cfam_elut <- merge(mer3_cfam_elut,mfurtt[,c("Protein.No.1","Transcript.No.1","Gene.No.1")],by.x="mfur.faa",by.y="Protein.No.1")
dim(mer4_cfam_elut) # not all of them have ferret orthologs

mer4_hsap_otters <- merge(mer3_hsap_otters,mfurtt[,c("Protein.No.1","Transcript.No.1","Gene.No.1")],by.x="mfur.faa",by.y="Protein.No.1")
dim(mer4_hsap_otters) # not all of them have ferret orthologs
mer4_cfam_otters <- merge(mer3_cfam_otters,mfurtt[,c("Protein.No.1","Transcript.No.1","Gene.No.1")],by.x="mfur.faa",by.y="Protein.No.1")
dim(mer4_cfam_otters) # not all of them have ferret orthologs

mer4_hsap_pbra <- merge(mer3_hsap_pbra,mfurtt[,c("Protein.No.1","Transcript.No.1","Gene.No.1")],by.x="mfur.faa",by.y="Protein.No.1")
dim(mer4_hsap_pbra) # not all of them have ferret orthologs
mer4_cfam_pbra <- merge(mer3_cfam_pbra,mfurtt[,c("Protein.No.1","Transcript.No.1","Gene.No.1")],by.x="mfur.faa",by.y="Protein.No.1")
dim(mer4_cfam_pbra) # not all of them have ferret orthologs
# now need to look up codeml ferret orthologs to see if they are in my VEP results

# elut:

vepgenesInCodeml_hsap_elut <- mer4_hsap_elut[mer4_hsap_elut$Gene.No.1 %in% allHighDomain_elut$Gene,] #50/250: 1227 ( was 204 (SNPs only); 1162 with indels)
dim(vepgenesInCodeml_hsap_elut)
vepgenesInCodeml_cfam_elut <- mer4_cfam_elut[mer4_cfam_elut$Gene.No.1 %in% allHighDomain_elut$Gene,] # (14; was 2 (snps only); 13 with indel ; so only gained 1)
dim(vepgenesInCodeml_cfam_elut)

# pbra: 
vepgenesInCodeml_hsap_pbra <- mer4_hsap_pbra[mer4_hsap_pbra$Gene.No.1 %in% allHighDomain_pbra$Gene,] 
dim(vepgenesInCodeml_hsap_pbra) # 1016 ( was 1012)
vepgenesInCodeml_cfam_pbra <- mer4_cfam_pbra[mer4_cfam_pbra$Gene.No.1 %in% allHighDomain_pbra$Gene,] 
dim(vepgenesInCodeml_cfam_pbra) # 13 ( was 13)
# cool, so elut gained genes (expected); pbra did not

## otters: need to see if they are pseudogenized in either elut or pbra (or both)
vepgenesInCodeml_hsap_otters_e <- mer4_hsap_elut[mer4_hsap_elut$Gene.No.1 %in% allHighDomain_elut$Gene,] # 
dim(vepgenesInCodeml_hsap_otters_e) # now 1227 (was 1162)
vepgenesInCodeml_hsap_otters_p <- mer4_hsap_pbra[mer4_hsap_pbra$Gene.No.1 %in% allHighDomain_pbra$Gene,] 
dim(vepgenesInCodeml_hsap_otters_p) # 1016 was 1012
vepgenesInCodeml_cfam_otters_e <- mer4_cfam_elut[mer4_cfam_elut$Gene.No.1 %in% allHighDomain_elut$Gene,] # 
dim(vepgenesInCodeml_cfam_otters_e) # 14 was 13
vepgenesInCodeml_cfam_otters_p <- mer4_cfam_pbra[mer4_cfam_pbra$Gene.No.1 %in% allHighDomain_pbra$Gene,] # 
dim(vepgenesInCodeml_cfam_otters_p) # 13 was 13

# note: the cfam / hsap separate lists are due to what the "transcript id" tag is for each cluster; most are from hsap, but some that are missing hsap have cfam as the TName. So that is why there are two lists.
############ FINAL WRITE OUT : ELUT #######
## Need to label these as "VEP-HIGH"
codemlElut$preVEPlabel <- codemlElut$label # save the old label for reference
codemlElut[codemlElut$TName %in% vepgenesInCodeml_cfam_elut$TName,]$label <- "VEP-HIGH"
codemlElut[codemlElut$TName %in% vepgenesInCodeml_hsap_elut$TName,]$label <- "VEP-HIGH"
sum(codemlElut$label=="VEP-HIGH") # 1241 ( was 1175)
# REMOVE vep-high genes: 
#write.table(codemlElut[codemlElut$label!="VEP-HIGH",],"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_elut/elut.RemovedVEPHighImpactGenes.putativelyReal.andUnexaminedGenes.20170922.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t")
# write them out separately:
write.table(codemlElut[codemlElut$label=="VEP-HIGH",],"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_elut/elut.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",row.names=F,quote=F,sep="\t")

############ FINAL WRITE OUT : PBRA #######
## Need to label these as "VEP-HIGH"
codemlPbra$preVEPlabel <- codemlPbra$label # save the old label for reference
codemlPbra[codemlPbra$TName %in% vepgenesInCodeml_cfam_pbra$TName,]$label <- "VEP-HIGH"
codemlPbra[codemlPbra$TName %in% vepgenesInCodeml_hsap_pbra$TName,]$label <- "VEP-HIGH"
sum(codemlPbra$label=="VEP-HIGH") # 1029 (was 1025)
# REMOVE vep-high genes: 
#write.table(codemlPbra[codemlPbra$label!="VEP-HIGH",],"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/pbra.RemovedVEPHighImpactGenes.putativelyReal.andUnexaminedGenes.20170922.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t")
# write them out separately:
write.table(codemlPbra[codemlPbra$label=="VEP-HIGH",],"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_pbra/pbra.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",row.names=F,quote=F,sep="\t")


############ FINAL WRITE OUT : OTTERS #######
## Need to label these as "VEP-HIGH"
# elut hsap/cfam Tnames:
codemlOtters$preVEPlabel <- codemlOtters$label # save the old label for reference

codemlOtters[codemlOtters$TName %in% vepgenesInCodeml_hsap_otters_e$TName,]$label <- "VEP-HIGH"
codemlOtters[codemlOtters$TName %in% vepgenesInCodeml_cfam_otters_e$TName,]$label <- "VEP-HIGH"

# pbra hsap/cfam Tnames:
codemlOtters[codemlOtters$TName %in% vepgenesInCodeml_hsap_otters_p$TName,]$label <- "VEP-HIGH"
codemlOtters[codemlOtters$TName %in% vepgenesInCodeml_cfam_otters_p$TName,]$label <- "VEP-HIGH"


sum(codemlOtters$label=="VEP-HIGH") # 1468 (was 1050 (elut) ; 1438 (elut + pbra) [note that some genes already missing from otters, so isn't a full 1631 that is the union of elut and pbra.  I checked the counts, looks ok]

# REMOVE vep-high genes: 
#write.table(codemlOtters[codemlOtters$label!="VEP-HIGH",],"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/otters.RemovedVEPHighImpactGenes.putativelyReal.andUnexaminedGenes.20170922.forPolysel.USETHIS.txt",row.names=F,quote=F,sep="\t")
# write them out separately:
write.table(codemlOtters[codemlOtters$label=="VEP-HIGH",],"/Users/annabelbeichman/Documents/UCLA/Otters/ComparativeGenomics/PamlMaterials/newCodemlResults_15spp/output_codeml_20170922/inspection_otters1/otters.VEPHighImpactGenes.50_250DPFilter.gotRemoved.txt",row.names=F,quote=F,sep="\t")


