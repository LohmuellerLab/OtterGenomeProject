
################ elut1 (gidget sso) #######################
elut_codingCallableSites_perChunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/calculatingCodingCallableSites/countCallableCDSRegions_50_250DPFilter/elut.callableSites.ferretCDSRegions.txt",header=T)
total_Coding_Callable_elut <- sum(elut_codingCallableSites_perChunk$callableCDSSites) 
total_Coding_Callable_elut #  30680418

################ pbra #######################

pbra_codingCallableSites_perChunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/calculatingCodingCallableSites/countCallableCDSRegions_50_250DPFilter/pbra.callableSites.ferretCDSRegions.txt",header=T)
total_Coding_Callable_pbra <- sum(pbra_codingCallableSites_perChunk$callableCDSSites) 
total_Coding_Callable_pbra # 30070013

############## elut2 (elfin NSO) #############
##### DON'T HAVE THIS FILE YET
elut2_codingCallableSites_perChunk <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/genome-stats/CallableCDSSites/elut2.callableSites.ferretCDSRegions.txt",header=T)
total_Coding_Callable_elut2 <- sum(elut2_codingCallableSites_perChunk$callableCDSSites) 
total_Coding_Callable_elut2 # 30432836

################# average ######################
average_Coding_Callable <- (total_Coding_Callable_elut+total_Coding_Callable_pbra+total_Coding_Callable_elut2)/3
average_Coding_Callable  # 30584207 
### write out a table:
allTotals <- rbind(total_Coding_Callable_elut,total_Coding_Callable_elut2,total_Coding_Callable_pbra,average_Coding_Callable)


write.table(allTotals,"/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/resultsCombining_all3Otters_inclNSO/VEP/TotalCodingCallableSites.and.averageCodingCallableSiteselut1pbra.50_250DPFilter.elut2.lib283.txt",row.names=T,quote=F,col.names=F)
