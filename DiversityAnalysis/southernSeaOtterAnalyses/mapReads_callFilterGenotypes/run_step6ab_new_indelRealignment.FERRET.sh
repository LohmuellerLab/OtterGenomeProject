#! /bin/bash
#$ -cwd
#$ -l h_rt=200:00:00,h_data=12G,highp,arch=intel*
#$ -o /u/scratch2/a/ab08028/otters/reports
#$ -e /u/scratch2/a/ab08028/otters/reports
#$ -m abe
#$ -M ab08028
# needs more than 50 hours!! check for completion.
source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
# usage: qsub -N indelRealign2steps <script> <samplePrefix: eg 01_Elut_CA_Gidget> 
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa 

# have to make these before 
IN_DIR=/u/scratch2/a/ab08028/otters/bams/
OUT_DIR=/u/scratch2/a/ab08028/otters/bams/IndelRealigned

PREFIX=$1


# step a: https://software.broadinstitute.org/gatk/gatkdocs/3.5-0/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php
# step b: Indel Realigner: https://software.broadinstitute.org/gatk/gatkdocs/3.7-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php
# note that they say you don't have to fix mates afterward anymore: http://gatkforums.broadinstitute.org/gatk/discussion/1562/need-to-run-a-step-with-fixmateinformation-after-realignment-step
# these filters are applied automatically:

    #MappingQualityZeroFilter
    #MalformedReadFilter
    #BadCigarFilter
    #BadMateFilter
    #UnmappedReadFilter
    #NotPrimaryAlignmentFilter
    #Platform454Filter
   	#FailsVendorQualityCheckFilter
   	#DuplicateReadFilter
    #MappingQualityUnavailableFilter
mkdir -p $OUT_DIR
# first step is to identify iffy regions with RealignerTargetCreator
java -Xmx8g -jar $GATK \
-T RealignerTargetCreator \
-R $REFERENCE \
-I $IN_DIR/${PREFIX}_Aligned_MarkDup.bam \
-o $OUT_DIR/${PREFIX}_Aligned_MarkDup_fromRTC.intervals

echo "done identifying regions, starting realignment now" 


# second step is to realign indels 

java -Xmx8g -jar $GATK \
-T IndelRealigner \
-R $REFERENCE \
-I $IN_DIR/${PREFIX}_Aligned_MarkDup.bam \
-targetIntervals $OUT_DIR/${PREFIX}_Aligned_MarkDup_fromRTC.intervals \
-o $OUT_DIR/${PREFIX}_Aligned_MarkDup_IndelRealigned.bam


