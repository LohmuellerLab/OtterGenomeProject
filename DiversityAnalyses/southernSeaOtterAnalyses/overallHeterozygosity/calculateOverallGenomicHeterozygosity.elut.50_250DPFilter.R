setwd("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/overallHeterozygosity")
# these results are from step_11_d (count vars)
# they are the counts of final HQ sites that have had all sites not passing filter REMOVED
# so no indels, no bad SNPs, all callable sites
# and then the final count of heterozygotes from my final variant files (after all filtering)
totalSites <- read.table("step-by-step-counts_50_250DPFilter/counts.9.HQSites.rmBadVars",header = F,stringsAsFactors = F)
# this is still wrong: 
totalHets <- read.table("step-by-step-counts_50_250DPFilter/counts.9.HQSites.rmBadVars.hets",header=F,stringsAsFactors = F)
totalHets2 <- read.table("step-by-step-counts_50_250DPFilter/counts.8.Vars.5.finalPASS.hets",header=F,stringsAsFactors = F) 
sum(totalHets2$V1) == sum(totalHets)

heterozygosity <- sum(totalHets$V1)/sum(totalSites$V1)
heterozygosity


cat(paste("Sea Otter Genome-Wide Heterozygosity: ",heterozygosity,"\n",sep=""),file="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/overallHeterozygosity/elut.GenomeWideHet.50_250DPFilter.txt")
### Also want to get total callable sites:

cat(paste("Genome_Wide_Callable ",sum(totalSites$V1),"\n","220_Scaffs_Callable ",sum(totalSites[1:220,"V1"]),sep=""),file="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/overallHeterozygosity/elut.TotalCallableSites.AllScaffs.220Scaffs.50_250DP.txt")
