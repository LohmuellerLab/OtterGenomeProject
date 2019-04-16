# submit two jobs: 
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts/
# note> there is no difference between PteBra_1.bam and PteBra_1_FastqToSam.bam 
# I just renamed it to be consistent with my pipeline
# they are just raw reads from fastq files

qsub -N step3 $scriptdir/run_step3_MarkIlluminaAdapters.AB.sh PteBra_1_FastqToSam.bam 


