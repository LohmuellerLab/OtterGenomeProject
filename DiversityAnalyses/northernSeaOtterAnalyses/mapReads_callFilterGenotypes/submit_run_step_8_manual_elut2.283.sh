scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/northern_sea_otter_mapping/TryWithJustOneLibrary_fromStep5On/IndelRealignment-BranchOfScripts
# usage:
# qsub script [pathtobamdir] [prefix] [L region number]
qsub $scriptdir/Step_8_a-alt_HaplotypeCaller-ScaffoldIntervals.lib283.sh /u/flashscratch/a/ab08028/otters/bams 02_Elut_SEAK_Elfin
# note the script is modified to use 02_Elut_SEAK_Elfin_Aligned_MarkDup_IndelRealigned.JustOneLib.alt.283_Filtered.bam
# as the input bam. this is in the regular ReadsFiltered directory, but just contains lib 283
# and has been indel-realigned. 
