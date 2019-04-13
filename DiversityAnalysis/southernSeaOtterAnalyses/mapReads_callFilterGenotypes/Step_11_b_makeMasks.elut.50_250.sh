#! /bin/bash
#$ -cwd
#$ -l h_rt=00:30:00,h_data=3G,highp
#$ -N makeMasks
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe
#$ -M ab08028
#$ -t 1-323
# pay attention to the diff between MASK (keeps what's in the mask) and negativemask (removes what's in the mask)
export PATH=$PATH:/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load python/3.4
module load perl
module load vcftools
module load samtools
module load bedtools

# date you called genotypes:
date=20171006
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter

i=${SGE_TASK_ID} # chunk
PREFIX=01_Elut_CA_Gidget
vcf=$wd/HQSitesOnly/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.vcf.gz # this has already gone through mask1 and negmask2

############## ONLY MAKE THE MASKS ONCE; one you run this script once, comment the making of the masks out.

########## Negative mask 1: variants that got filtered out 
echo "making negative mask 1 (variants that got filtered out)"
mkdir -p $wd/variants/badVarsNegativeMasksForMSMC/

vcf-isec -c $wd/variants/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.1.Variants.vcf.gz $wd/variants/${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.vcf.gz \
| grep -v "#" | awk '{OFS="\t";print $1,$2-1,$2}' > $wd/variants/badVarsNegativeMasksForMSMC/chunk_${i}_badVars.0based.bed 

# merge it
bedtools merge -i $wd/variants/badVarsNegativeMasksForMSMC/chunk_${i}_badVars.0based.bed > $wd/variants/badVarsNegativeMasksForMSMC/chunk_${i}_badVars.0based.merged.bed


##### Negative mask 3: sites with multinucleotide reference allele (so weird)
echo "making negative mask 3 (note we skipped mask2, HQ sites already filtered) -- sites where reference allele is multinucleotide"
mkdir -p $wd/multiNucRefAlleleBedCoords
zcat $vcf | grep -v "#" | awk '{OFS="\t";if(length($4)>1)print $1,$2-1,$2}' > $wd/multiNucRefAlleleBedCoords/chunk_${i}.multinuc.RefAllele.0based.bed

# merge it:
bedtools merge -i $wd/multiNucRefAlleleBedCoords/chunk_${i}.multinuc.RefAllele.0based.bed > $wd/multiNucRefAlleleBedCoords/chunk_${i}.multinuc.RefAllele.0based.merged.bed

# regions you want to EXCLUDE (indels, any snps that didn't pass filters) 
sleep 5m
