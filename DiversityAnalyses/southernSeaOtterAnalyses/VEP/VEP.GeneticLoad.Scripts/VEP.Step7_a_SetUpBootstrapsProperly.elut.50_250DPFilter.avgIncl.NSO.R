######### This scripts sets up a file with the numbers of derived alleles found in the NS, S and SG states in your data 
spp="elut"
prefix="01_Elut_CA_Gidget"
date="20171006"
datadir="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/VEP/VEP-50_250DPFilter_canon_cds/includesShared/" ### make specific to elut or pbra 
dir.create(paste(datadir,spp,"_derivedAllele_bootstraps/",sep=""))
inputFile=paste(datadir,prefix,".",date,".VepResultsCountSummary.txt",sep="")
# this file is the same for all 3 spp:
########## UPDATED 20190201 to use the average from all three otter species:  ############

# *** 20180201: this is a new averages file (product of getAverages /Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mapping_scripts_50_250DPFILTER_elutpbra_USETHIS/CallableCDSSites/getAverageCodingCallableSites.forRescalingLoadCalculations.elut1.elut2.pbra.R): *********
totalsCodingCallable <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/resultsCombining_all3Otters_inclNSO/VEP/TotalCodingCallableSites.and.averageCodingCallableSiteselut1pbra.50_250DPFilter.elut2.lib283.txt",header=F)
average_Coding_Callable <- totalsCodingCallable[totalsCodingCallable$V1=="average_Coding_Callable",]$V2
nsites=average_Coding_Callable # number of sites to draw out for each bootstrap (with replacement, use average callables sites in coding regions)
# new average for all 3 otters 30394422
# elut and pbra totals:
elut_total_coding_callable = totalsCodingCallable[totalsCodingCallable$V1=="total_Coding_Callable_elut",]$V2
pbra_total_coding_callable = totalsCodingCallable[totalsCodingCallable$V1=="total_Coding_Callable_pbra",]$V2
# added elut2:
elut2_total_coding_callable = totalsCodingCallable[totalsCodingCallable$V1=="total_Coding_Callable_elut2",]$V2
### get the appropriate total for your species: 
if(spp=="pbra"){
  totalCoding=pbra_total_coding_callable
} else if(spp=="elut"){
  totalCoding=elut_total_coding_callable
} else if(spp=="elut2"){
  totalCoding=elut2_total_coding_callable
} else {
  print("something has gone wrong! check your species!")
}
# check: 
totalCoding==elut_total_coding_callable

##################### TOTALS ###########
################### overall hets/homs category ######################
# these are amounts of overal homsalt/hets/homsref in coding sites 
############# pbra updated to be 50/250 DP filter
########## THESE ARE FOR ELUT: 
####### These are totals of all coding sites (Hom alt , hom ref, hets), regardless of VEP annotation
# some are unannotated by VEP 
coding_homAlt_perchunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/calculatingCodingCallableSites/countCallableCDSRegions_50_250DPFilter/elut.callableSites.ferretCDSRegions.VariantHoms.txt",header=T)
coding_homAlt_total <- sum(coding_homAlt_perchunk$callableCDSSites)
coding_homAlt_total

coding_homRef_perchunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/calculatingCodingCallableSites/countCallableCDSRegions_50_250DPFilter/elut.callableSites.ferretCDSRegions.NonVariant.HomRefs.txt",header=T)
coding_homRef_total <- sum(coding_homRef_perchunk$callableCDSSites)
coding_homRef_total

coding_het_perchunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/calculatingCodingCallableSites/countCallableCDSRegions_50_250DPFilter/elut.callableSites.ferretCDSRegions.VariantHets.txt",header=T)
coding_het_total <- sum(as.numeric(coding_het_perchunk$callableCDSSites))
coding_het_total

# check totals add up:
elut_total_coding_callable == coding_het_total + coding_homAlt_total + coding_homRef_total # TRUE 
####### Need to create your input dataset of all sites with correct numbers of derived alleles: 
####### making a dummy set of sites that matches the counts of all the different categories of sites you annotated (synonymous (hom/het), missense (hom/het), SG (hom/het) and unnanotated)
# (it's basically just your dataset, but highly simplified)
# going to resample the avg number of called sites (across the 3 otters) form this dataset with replacement to make your bootstraps
## Use count summaries to make the input file:
# for testing:
inputCounts <- read.table(inputFile,header=T)
vep_Mis_het <- inputCounts[inputCounts$category=="missense" & inputCounts$genotype=="Heterozygous",]$count # missense het count
# make a dataframe of sites with that count and category (derived allele count is 1 for hets)
df_Mis_het <- data.frame(Annot=rep("missense",vep_Mis_het),DerivedCount=rep(1,vep_Mis_het))

vep_Mis_hom <- inputCounts[inputCounts$category=="missense" & inputCounts$genotype=="HomozygousAlt",]$count
# derived allele count is 2 for homs
df_Mis_hom <- data.frame(Annot=rep("missense",vep_Mis_hom),DerivedCount=rep(2,vep_Mis_hom))

vep_Syn_het <- inputCounts[inputCounts$category=="synonymous" & inputCounts$genotype=="Heterozygous",]$count
df_syn_het <- data.frame(Annot=rep("synonymous",vep_Syn_het),DerivedCount=rep(1,vep_Syn_het))

vep_Syn_hom <- inputCounts[inputCounts$category=="synonymous" & inputCounts$genotype=="HomozygousAlt",]$count
df_syn_hom <- data.frame(Annot=rep("synonymous",vep_Syn_hom),DerivedCount=rep(2,vep_Syn_hom))

vep_SG_het <- inputCounts[inputCounts$category=="stopgained" & inputCounts$genotype=="Heterozygous",]$count
df_sg_het <- data.frame(Annot=rep("stopgained",vep_SG_het),DerivedCount=rep(1,vep_SG_het))

vep_SG_hom <- inputCounts[inputCounts$category=="stopgained" & inputCounts$genotype=="HomozygousAlt",]$count
df_sg_hom <- data.frame(Annot=rep("stopgained",vep_SG_hom),DerivedCount=rep(2,vep_SG_hom))

############# Unannotated bases check this carefully ############
#totalUnannotated <- totalSiteCount - totalAnnotated
totalUnannotedHet <- coding_het_total - (vep_Mis_het+vep_Syn_het+vep_SG_het)

totalUnannotatedHomAlt <- coding_homAlt_total - (vep_Mis_hom+vep_Syn_hom+vep_SG_hom)

df_totalUnannotedHet <- data.frame(Annot=rep("Not_Annotated",totalUnannotedHet),DerivedCount=rep(1,totalUnannotedHet))

df_totalUnannotedHomAlt <- data.frame(Annot=rep("Not_Annotated",totalUnannotatedHomAlt),DerivedCount=rep(2,totalUnannotatedHomAlt))
############### Homozygous reference ############
df_HomRef <- data.frame(Annot=rep("HomozygousReference",coding_homRef_total),DerivedCount=rep(0,coding_homRef_total)) # this should be derived allele count of 0 because it is homref!
##### Make input dataframe: #######
# all categories from vep, all unannotated, all hom ref (should equal 100K?)
inputData <- rbind(df_Mis_het,df_Mis_hom,df_syn_het,df_syn_hom,df_sg_het,df_sg_hom,df_totalUnannotedHet,df_totalUnannotedHomAlt,df_HomRef)
dim(inputData)[1]==elut_total_coding_callable # should equal total coding sites. Good!
# write this out:
write.table(inputData,paste(datadir,spp,"_derivedAllele_bootstraps/",spp,".forBootstraps.StrippedDownDataset.AllCallableSites.derivedCountfromVEP.Annotation.50_250DPFilter.AvgIncl.NSO.283.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
### now do this for elut. then run boots
