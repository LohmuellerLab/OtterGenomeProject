setwd("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/genome-stats/overall-heterozygosity")
# these results are from step_11_d (count vars)
# they are the counts of final HQ sites that have had all sites not passing filter REMOVED
# so no indels, no bad SNPs, all callable sites
# and then the final count of heterozygotes from my final variant files (after all filtering)

totalSites <- read.table("step-by-step-counts/counts.9.HQSites.rmBadVars",header = F)
totalHets <- read.table("step-by-step-counts/counts.9.HQSites.rmBadVars.hets",header=F) 
# use this; most reliable in case masking went weird. 
totalHets2 <- read.table("step-by-step-counts/counts.8.Vars.5.finalPASS.hets",header=F) 
sum(totalHets$V1 == totalHets2$V1) # 323 good. means dimensions are the same
sum(totalHets) == sum(totalHets2) # TRUE means that count of hets from two different files are the same. good.
heterozygosity <- sum(totalHets$V1)/sum(totalSites$V1)
# 0.00030824 northern sea otter
heterozygosity
cat(paste("Northern Sea Otter (elut2) Genome-Wide Heterozygosity: ",heterozygosity,"\n",sep=""),file="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/genome-stats/overall-heterozygosity/elut2.GenomeWideHet.50_250DPFilter.283.txt")



### Also want to get total callable sites:
cat(paste("Genome_Wide_Callable ",sum(totalSites$V1),"\n","220_Scaffs_Callable ",sum(totalSites[1:220,"V1"]),sep=""),file="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/genome-stats/overall-heterozygosity/elut2.TotalCallableSites.AllScaffs.220Scaffs.50_250DP.283.txt")
