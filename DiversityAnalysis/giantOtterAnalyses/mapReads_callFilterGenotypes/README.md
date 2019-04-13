# Guide to read mapping, genotype calling and filtering  
*based on GATK best practices (Van der Auwera et al. 2013)*  
*full details in SI Methods of Beichman et al (2019)*  
- **Step 0:** download data  
- **Step 1:** Quality assessment: FASTQC (Andrews 2010) was used to assess the quality of sequencing reads (script not included)  
- **Step 2:** Fastq to uBAM: Fastq files were converted to uBAM files using Picard FastqToSam  
- **Step 3:** Mark adapters: Illumina adapters were marked using Picard MarkIlluminaAdapters and clipped using Picard FastqToSam  
- **Step 4:** Align reads 

  **a.**	Clip adapters: Illumina adapters were converted to Ns using Picard FastqToSam  
  
  **b.**	Reads were aligned to the reference genome (MutPutFur1.0) using BWA-MEM  
  
  **c.**	Information from the unmapped bam file is added to the alignment using Picard MergeBamAlignment, with the MostDistant alignment strategy, with any amount of insertions or deletions allowed.  
  
- **Step 5: Mark duplicates: duplicate reads were marked using Picard MarkDuplicates  
- **Step 6: Indel realignment: indels were realigned using GATKÕs RealignerTargetCreator and IndelRealigner  
a+b: indel realignment (two steps; first takes ~6h, second takes a *long time* (3-7 days!)) this is a legacy step to match with past data. Hap Caller *does* indel realignment, so this is not needed. However, I did it for my first two otters, so I did it with any data I wanted to compare with them, just to be safe.  
- **Step 7**: Removal of bad reads: using SAMtools view (Li et al. 2009), reads are removed that are not properly aligned (-f 2), are secondary alignments (-F 256), have a mapping quality score of less than 30, or are PCR or optical duplicates (-F 1024).  
- **Step 8**: Generate gVCFs: variants were called using GATKÕs HaplotypeCaller to generate a gVCF file for each of the largest 224 scaffolds separately, and the remaining ~7,500 smaller scaffolds in groups of 76 (to allow for parallelization). All sites passing a base quality score of 20 were called, whether or not they were variant.  
- **Step 9**: Call genotypes: genotypes were called using GATKÕs GenotypeGVCFs to generate vcf files.  
- **Step 10**: Get callable sites:  

  **a-d.** get the distribution of read depth per site (DP) was calculated from the largest 50 scaffolds, calculate mean DP  
  
  **e-f.** do the filtering based on DP filter  
  
- **Step 11:** Variant filtering  

  **a-c:**	Variants sites were selected from the overall vcf file and were filtered according to GATK recommendations (Van der Auwera et al. 2013): filtering out variants with quality-by-depth score (QD) < 2.0, with Fisher Strand bias (FS) > 60, mapping quality (MQ) < 40, z-score of rank sum test of alternate vs. reference read mapping qualities (MQRankSum) < -12.5, or z-score of alternate vs. reference read position bias (ReadPosRankSum) less than -8.0.   
Biallelic SNPs were selected, and SNPs clustered with more than 3 SNPs in 10bp were removed.  

  **d:** count variants at each stage of filtering  

- **Step 12:** Indel filtering  
