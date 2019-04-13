## DiversityAnalyses
### *Scripts to analyze genetic diversity, demography and genetic load*

#### giantOtterAnalyses / southernSeaOtterAnalyses / northernSeaOtterAnalyses:
*Population genetics analyses for each taxa*
- **mapReads_callFilterGenotypes**: scripts
to carry out QC, align sequencing reads, and filter genotypes.
Modified from pipelines developed by Dr. Jacqueline Robinson, 
Dr. Clare Marsden and Dr. Tanya Phung. Has its own, more detailed, readme in the directory.
- **overallHeterozygosity**: calculate genome-wide heterozygosity
- **slidingWindowHeterozygosity**: calculate heterozygosity in sliding windows 
(modified from Dr. Jacqueline Robinson)
- **CallableCDSSites**: count callable sites in coding regions of the ferret genome (cds)
- **VEP**: run Ensembl's variant effect predictor (McLaren et al. 2016). 
Use it to detect potentially pseudogenized genes (VEP.Pseuodgenes.Scripts)
and synonymous / missense / stop-gained variants (VEP.GeneticLoad.Scripts).  
The genetic load scripts run the variant effect predictor on your sites, 
filters the output by synonymous/missense/stop-gained 
and helps you set up bootstraps to compare burden of 'deleterious' mutations . 
This analysis is highly sensitive to problems in the data, so it is important that levels 
of synonymous mutations match across taxa.
- **ROH**: detect runs of homozygosity 
- **MSMC_DemographicInference**: carry out MSMC inference and run simulations
to check whether the inferred demography predicts heterozygosity

#### analysesComparingSpecies:
*Comparisons across southern/northern sea otter and giant otter*
- **combineVCFs**: combine genotype vcf files across all three taxa
- **averageCodingCallableSites**: calculate average called coding sites across taxa
- **ABBA-BABA**: check species topology and detect admixture (or lack thereof) between taxa
- **compareVEPResults_betweenElutPbra**: compare southern sea otter and giant otter VEP sites
- **mutationRate**: calculate average mutation rate between southern sea otter and giant otter
