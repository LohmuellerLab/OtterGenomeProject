# 5c: add outgroups (Niimura), human sequence, and Gang Li's representative sequences together with your result
spp=$1
dir=/work2/abeichman/gangLi_OlfactionPipeline_${spp}
cat $dir/Human_OR2J3.fasta \
$dir/niimura.genbank.outgroupSeqs.fasta \
$dir/olfac_5sp_2seq_query.fasta \
$dir/step_5_chooseM/step_5_result.pickedM.fasta > \
$dir/step_5_chooseM/step_5_result.pickedM.wHuman.wOutgroups.wRepresentatives.fasta
