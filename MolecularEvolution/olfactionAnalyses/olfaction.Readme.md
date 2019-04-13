## Guide to olfactory gene detection:
#### *based on scripts from Dr. Gang Li and protocols from Niimura (2013) and Montague et al. (2014)*
#### Details in SI methods of Beichman et al. (2019)

#### representative sequences used to find ORGs: 
- **olfac_5sp_2seq_query.fasta**: representative sequences of olfactory receptors from Montague et al. (2014)
- **Human_OR2J3.fasta**: human representative ORG sequence
- **niimura.genbank.outgroupSeqs.fasta**: outgroup sequences (non-olfactory receptor genes) from Niimura (2013)

#### FunctionalOlfactoryReceptorGene_Analysis
*Detecting functional ORGs*

**ferret_ORGs / giant_otter_ORGs / northern_sea_otter_ORGs / southern_sea_otter_ORGs:**

*scripts to detect functional ORGs in each taxa's genome*
- **step_1_ab_blastDb_blastResults**: generate a blast database and tblastn the sequences from 5 mammal species against the genome. 
- **step_1_c_CleanedByEvalue_ABScript**: clean blast results and select best matches with a length filter: removes hits <250 AAs, removes overlapping hits by picking the one with the best evalue
- **step_2_getSequencesFromGenome**: pull the blast results from the genome +- 300bp on either side 
- **step_3_findORFs**: detect open reading frames (run twice, once with and once without print statements to get full info on alignments)
- **step_4_alignWithMafft_inclHuman**: align detected ORG sequences, visually inspect alignments to remove sequences failing Niimura's critera for functional ORGs
- **step_5_chooseM**: detect coordinates of start-codon
- **step_6_alignWithOutgroups**: align to outgroup non-ORG sequences
- **step_7_plotinROnDesktop**: plot in phylogenetic tree
- **step_8_alignClassesSeparately**: detect and align Class I and Class II ORGs separately
- **step_9_finalFiltering**: do last checks for any further ORGs failing requirements from Niimura (2013)
- **step_10_getFinalGenomeCoordsAndSeqs**: get coordinates of genome of final set of functional ORGs and their sequences in fasta format

#### PseudogenizedOlfactoryReceptorGenes_Analysis
*Detecting pseudogenized ORGs*

**ferret_pseudo_ORGs / giant_otter_pseudo_ORGs / northern_sea_otters_pseudo_ORGs / southern_sea_otter_pseudo_ORGs:**

*scripts to detect pseudogenized ORGs in each taxa's genome*

- **step 1**: clean up blast results without imposing a length restriction (because want to detect fragments)
- **step 2**: get fasta sequences of blast results 
- **step 3**: find open reading frames
- **step 4**: align sequences and visually categorize into whether they are missing N and/or C terminal ends
- **step 5**: process and label pseudogenes in R
