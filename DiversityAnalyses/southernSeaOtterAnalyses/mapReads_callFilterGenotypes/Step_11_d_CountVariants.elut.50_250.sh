#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=5G,highp
#$ -N countSites
#$ -o /u/flashscratch/a/ab08028/otters/reports/
#$ -e /u/flashscratch/a/ab08028/otters/reports/
#$ -m abe

date=20171006
wd=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter
out=$wd/step-by-step-counts
mkdir -p $out
### Where am I losing SNPs: the script:
> $out/counts.1.raw.genotypes
> $out/counts.2.HQ.Sites.wDotGTs
> $out/counts.3.HQ.Sites.rmDotGTs
> $out/counts.4.Vars.1.All
> $out/counts.5.Vars.2.HF
> $out/counts.6.Vars.3.BiSNP
> $out/counts.7.Vars.4.rmClust
> $out/counts.8.Vars.5.finalPASS
> $out/counts.9.HQSites.rmBadVars
> $out/counts.1.raw.genotypes.hets
> $out/counts.2.HQ.Sites.wDotGTs.hets
> $out/counts.3.HQ.Sites.rmDotGTs.hets
> $out/counts.4.Vars.1.All.hets
> $out/counts.5.Vars.2.HF.hets
> $out/counts.6.Vars.3.BiSNP.hets
> $out/counts.7.Vars.4.rmClust.hets
> $out/counts.8.Vars.5.finalPASS.hets
> $out/counts.8.Vars.5.finalPASS.homs
> $out/counts.9.HQSites.rmBadVars.hets

for i in {1..323}
do
echo $i
# 1. Raw genotypes:
zcat $wd/raw_variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.vcf.gz | grep -c -v "#"  >> $out/counts.1.raw.genotypes
zcat $wd/raw_variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.vcf.gz | grep -c "0/1"  >> $out/counts.1.raw.genotypes.hets

# 2. HQ Sites (within DP range), but didn't remove dotted genotypes
zcat $wd/HQSitesOnly-stillbadgenotypes/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.vcf.gz | grep -c -v "#" >> $out/counts.2.HQ.Sites.wDotGTs
zcat $wd/HQSitesOnly-stillbadgenotypes/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.vcf.gz | grep -c "0/1" >> $out/counts.2.HQ.Sites.wDotGTs.hets

# 3. HQ Sites, remove dotted GTs
zcat $wd/HQSitesOnly/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.rmDotGenotypes.vcf.gz | grep -c -v "#" >> $out/counts.3.HQ.Sites.rmDotGTs
zcat $wd/HQSitesOnly/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.rmDotGenotypes.vcf.gz | grep -c "0/1"  >> $out/counts.3.HQ.Sites.rmDotGTs.hets

# 4. Pull out all variant sites (1/1, 0/1, indels and snps)
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.1.Variants.vcf.gz | grep -c -v "#" >> $out/counts.4.Vars.1.All
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.1.Variants.vcf.gz | grep -c "0/1" >> $out/counts.4.Vars.1.All.hets

# 5. Hard filters
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.2.HF.Variants.vcf.gz | grep -v "#" | grep -c PASS >> $out/counts.5.Vars.2.HF
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.2.HF.Variants.vcf.gz | grep -v "#" |grep PASS | grep -c "0/1" >> $out/counts.5.Vars.2.HF.hets

# 6. Biallelic SNPs: 
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.3.BiSNP.HF.Variants.vcf.gz | grep -v "#" | grep -c PASS >> $out/counts.6.Vars.3.BiSNP
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.3.BiSNP.HF.Variants.vcf.gz | grep -v "#" | grep PASS | grep -c "0/1" >> $out/counts.6.Vars.3.BiSNP.hets

# 7. remove clusters:
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.4.rmClust.BiSNP.HF.Variants.vcf.gz | grep -v "#" | grep -c PASS >> $out/counts.7.Vars.4.rmClust
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.4.rmClust.BiSNP.HF.Variants.vcf.gz | grep -v "#" | grep PASS | grep -c "0/1" >> $out/counts.7.Vars.4.rmClust.hets 

# 8. final vars
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.vcf.gz | grep -v "#" | grep -c PASS >> $out/counts.8.Vars.5.finalPASS
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.vcf.gz | grep -v "#" | grep PASS | grep -c "0/1" >> $out/counts.8.Vars.5.finalPASS.hets
zcat $wd/variants/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.5.PASS-ONLY.rmClust.BiSNP.HF.Variants.vcf.gz | grep -v "#" | grep PASS | grep -c "1/1" >> $out/counts.8.Vars.5.finalPASS.homs

# 9. count final HQ sites after masking bad variants:
zcat $wd/HQSitesOnly-MaskedBadVariantsIndels/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz | grep -c -v "#" >> $out/counts.9.HQSites.rmBadVars
zcat $wd/HQSitesOnly-MaskedBadVariantsIndels/01_Elut_CA_Gidget.raw_variants.${i}.20171006.HQsites.Only.rmDotGenotypes.rmBadVars.vcf.gz | grep -c "0/1" >> $out/counts.9.HQSites.rmBadVars.hets

# so genome wide het will be 9a/9b
done
