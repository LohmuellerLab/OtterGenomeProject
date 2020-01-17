# have to manually replace the human sequence (not output for some reason...)
input=step_5_result_mafftAlignment.pickedM.fasta
human=humanSeqOriginal.fasta
outgroups=niimura.genbank.outgroupSeqs.fasta
grep -v ">Human_OR2J3" $input > ${input%.fasta}.noHuman.fasta
cat $human $outgroups ${input%.fasta}.noHuman.fasta > ${input%.fasta}.wHuman.wOutgroups.fasta

