# Sea otter and giant otter genome project scripts
*These scripts were used to generate the analyses for Beichman et al. (2019).*  
*Details of the methods are in the supplementary materials of the paper.*  
*Contact: Annabel Beichman <annabel.beichman[at]gmail.com>*  
Note: these scripts are internal to this project, 
based on my file systems and directory organization.
Feel free to adapt and reuse in your own work with the caveat 
that they are not currently designed for generic use,
so results are not guaranteed.

## MolecularEvolution: 
### *Scripts to analyze positive selection and gene loss*
#### proteinOrtho: 
*Scripts to find orthologous protein sequences across mammal species
using ProteinOrtho (Lechner et al. 2011)*

#### codemlAnalyses: 
*Scripts to carry out multi-species sequence alignments, 
filtering and masking of alignments, 
and tests for positive selection 
using the codeml branch-site test (Yang et al. 2007).
Has its own, more detailed, readme in the directory.*

#### polyselAnalyses:
*Scripts to carry out polygenic selection analysis 
on results from codeml using scripts and guides 
from Daub et al. (2016) https://github.com/CMPG/polysel*

#### pseudoPipeAnalyses:
*Scripts to detect pseudogenes using PseudoPipe (Zhang et al. 2006)*

#### olfactionAnalyses:
*Scripts to detect functional and pseudogenized 
olfactory receptor genes (ORGs), 
scripts modified from Dr. Gang Li.
Has its own, more detailed, readme in the directory.*

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

## PlottingScripts:
### *Scripts to generate several plots from Beichman et al. (2019)*
- **HetIUCNComparison**: plot genome-wide heterozygosity based on values from the literature (modified from Robinson et al, 2016)
- **ROH**: plot runs of homozygosity
- **SlidingWindowHet**: plot sliding window heterozygosity
- **MSMC**: plot MSMC trajectories and simulations


