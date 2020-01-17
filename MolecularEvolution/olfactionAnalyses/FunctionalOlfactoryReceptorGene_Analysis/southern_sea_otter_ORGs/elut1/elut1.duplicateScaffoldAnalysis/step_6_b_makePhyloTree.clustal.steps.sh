## command line clustal tree

clustalw=/home/abeichman/bin/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2
# says to use kimura's correction, remove sites with gaps, 1000 bootstraps, and bootstrap labels on node not branch

$clustalw \
-CLUSTERING=NJ -bootstrap=1000 \
-KIMURA -TOSSGAPS -BOOTLABELS=node \
-infile=step_6_result_mafftAligment.wOutgroups.20171103.fasta \
-outfile=step_6_result_mafftAligment.wOutgroups.phb
