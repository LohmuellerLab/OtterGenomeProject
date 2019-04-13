# Want to get CDS count per gene
# start with cds per 

module load bedtools
# these are all cds genotypes (need total not cdsVCF
date=20171006
vcfdir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_${date}_50_250DPFilter/cds_vcfs/
PREFIX=01_Elut_CA_Gidget
# get a header:
grep "#" ${vcfdir}/${PREFIX}.raw_variants.1.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf > ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.AllCDSSites.vcf

for i in {1..323}
do
echo $i
vcf=${PREFIX}.raw_variants.${i}.${date}.HQsites.Only.rmDotGenotypes.rmBadVars.CDSOnly.fromMergedBed.vcf
# get ALL sites (including hom ref):
grep -v "#" ${vcfdir}/${vcf} >> ${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.AllCDSSites.vcf
done

allCDSSites=${vcfdir}/${PREFIX}.raw_variants.cdsSequence.AllScaffsConcat.${date}.HQsites.Only.AllCDSSites.vcf
gff=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/MusPutFuro1.0.91.cdsOnly.gff # this is the ferret gff file
# this is a gff file of all CDS sequences in the ferret genome 
mkdir -p $vcfdir/getCDSCountPerGene
bedtools intersect -c -a $gff -b $allCDSSites > $vcfdir/getCDSCountPerGene/${PREFIX}.${date}.cds.countsPerExon.txt
### want the vcf to be all sites 
# so the goal of this is to count the number of hits in the site by site VCF (that includes ALL cds sites, including hom ref)
# that are part of each entry in A.
# this will be exon by exon; so then need to use dplyr to get total per gene (hopefully works)

#  	For each entry in A, report the number of hits in B while restricting to -f. Reports 0 for A entries that have no overlap with B. Restricted by -f and -r.
