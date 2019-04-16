# submit two jobs: 
# gidget lanes 2 and 6
qsub -N step3.1 /u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts/run_step3_MarkIlluminaAdapters.AB.sh RWAB001_L002_FastqToSam.bam


qsub -N step3.2 /u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts/run_step3_MarkIlluminaAdapters.AB.sh RWAB001_S6_L006_FastqToSam.bam
