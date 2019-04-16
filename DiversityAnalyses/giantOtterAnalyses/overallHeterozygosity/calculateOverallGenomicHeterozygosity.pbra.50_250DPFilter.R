setwd("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/overall-heterozygosity")
# these results are from step_11_d (count vars)
# they are the counts of final HQ sites that have had all sites not passing filter REMOVED
# so no indels, no bad SNPs, all callable sites
# and then the final count of heterozygotes from my final variant files (after all filtering)
# 20171208: these are temp files because only go up to ~224; need for all scaffolds (waiting on 225,226,227)  then redo counts
totalSites <- read.table("step-by-step-counts_50_250DPFilter/counts.9.HQSites.rmBadVars",header = F)
totalHets <- read.table("step-by-step-counts_50_250DPFilter/counts.9.HQSites.rmBadVars.hets",header=F) 
# use this; most reliable in case masking went weird. 
totalHets2 <- read.table("step-by-step-counts_50_250DPFilter/counts.8.Vars.5.finalPASS.hets",header=F) 
sum(totalHets == totalHets2) # 323 good. 
sum(totalHets) == sum(totalHets2) # TRUE
heterozygosity <- sum(totalHets$V1)/sum(totalSites$V1)
# 0.0006282961 with 50/250 dp filter 
heterozygosity
cat(paste("Giant Otter Genome-Wide Heterozygosity: ",heterozygosity,"\n",sep=""),file="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/overall-heterozygosity/pbra.GenomeWideHet.50_250DPFilter.txt")



### Also want to get total callable sites:
cat(paste("Genome_Wide_Callable ",sum(totalSites$V1),"\n","220_Scaffs_Callable ",sum(totalSites[1:220,"V1"]),sep=""),file="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/overall-heterozygosity/pbra.TotalCallableSites.AllScaffs.220Scaffs.50_250DP.txt")
