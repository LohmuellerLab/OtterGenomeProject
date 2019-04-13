scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/northern_sea_otter_mapping/TryWithJustOneLibrary_fromStep5On



# I am giving output name : arbitrary Number_Population_ID_

qsub -N MarkDup $scriptdir/run_step5_MarkDuplicates.AB.elut2.JustOneLib.sh SRX2967283_Aligned.bam 02_Elut_SEAK_Elfin_Aligned_MarkDup.JustOneLib.alt.283.bam

