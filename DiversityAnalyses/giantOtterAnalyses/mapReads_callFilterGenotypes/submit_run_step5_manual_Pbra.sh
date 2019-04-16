SCRIPT_DIR=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts

# I am giving output name : arbitrary Number_Population_ID_

qsub -N MarkDup $SCRIPT_DIR/run_step5_MarkDuplicates.Pbra.AB.sh PteBra_1_Aligned.bam PteBra_1_Aligned_MarkDup.bam

