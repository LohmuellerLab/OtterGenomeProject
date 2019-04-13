# Guide to codeml analysis  

*Scripts to align 1:1 orthologs, filter the alignments, and detect selection using codeml's branch-site test* 
*Details in SI methods of Beichman et al. (2019)*

- **step 1**: evaluate results from proteinOrtho and pull out 1:1 orthologs  
- **step 2**: set up directory structure for codeml analyses  
- **step 3**: get transcript sequences for each cluster of 1:1 orthologs  
- **step 4**: alignment  
  **a.** run GUIDANCE2 with PRANK to align sequences and detect low-quality sequences  
  **b.** remove bad sequences (non-div by 3, premature stop codon, low quality) and re-run Guidance (iterate this step until no more failures)  
  **c.** Mask codons flagged by GUIDANCE2 as unreliably aligned, convert alignment to phylip format  
- **step 5**: filter alignments:  
  **a.** prune tree topologies to remove species whose sequences were removed in step 4  
  **b.** CHOOSE A FILTERING METHOD (or skip): Gblocks block based filtering or SWAMP sliding window filtering  
- **step 6**: run codeml's branch-site test (three replicates)  
- **step 7**: Gather codeml results together into a table  
- **step 8**: Analyse and summarise codeml results  
  **a-i.** prepare information about all the alignment clusters  
  **a-ii.** process results, do p-value correction for multiple testing across all branch-genes and make list of genes past significance threshold  
  **b.** gather up significant genes so that you can visually inspect their alignments.  
  *< visually inspect all BEB sites within significant alignments, flag alignments that are suspect and remove them from dataset >**  
  *< intersect with results from VEP based on reads mapped to ferret genome; genes with high-impact VEP sites should also be removed as they may be pseudogenized >**   
  **c.** after removing all suspect alignments, reprocess data  
  
 