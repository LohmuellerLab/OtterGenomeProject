SCRIPT_DIR=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/scripts

# I am giving output name : arbitrary Number_Population_ID_

qsub -N MarkDup $SCRIPT_DIR/run_step5_MarkDuplicates.AB.sh RWAB001_L002_Aligned.bam  RWAB001_S6_L006_Aligned.bam 01_Elut_CA_Gidget_Aligned_MarkDup.bam

