#! /bin/bash
#$ -l h_rt=01:00:00,h_data=1G
#$ -cwd
#$ -N get_seqs
#$ -m abe
#$ -M ab08028
#$ -o prankReports
#$ -e prankReports
#$ -t 1-77

seqtk=/u/home/a/ab08028/klohmueldata/annabel_data/bin/seqtk/seqtk
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file 
allFasta=all15spp_CDS.TranscriptIDs.1.20170606.fna # make a combined fasta file of all of each species cds sequences of genes with associated transcript IDs; make sure there are no UTRs.


j=`printf "%02d" $((10#$SGE_TASK_ID - 10#1))`
start=$((10#$j*200))
echo "start: " $start
end=$((10#$j*200+199))
echo "end: " $end
for (( i=$start; i<=$end; i++ ))
do 
echo $i
# for awk need i plus 1 (zeroth cluster is in first column)
p=$(($i+1))
# get human T ID to be overall gene name
TName=`awk -v a=${row1} 'FNR==a {print $"'$p'"}' $clusters`
# but if TName is missing in human, sub in the dog:
if [[ $TName = "NA" ]]
then 
TName=`awk -v a=${row2} 'FNR==a {print $"'$p'"}' $clusters`
fi

cluster=`printf "%05d" $i`

# seqtk works well : if a sequence in the list cannot be found, it just doesn't add it
# so anything called NA or UMAR_NA or OROS_NA will just not be added to the final transcript file

$seqtk subseq $allFasta group_${j}/cluster_${cluster}_group_${j}_${TName}/TrIDS_cluster_${cluster}_group_${j}_${TName}.txt > group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna

# then want to replace name of transcript with the species name (-i'' means to do it in place)
# amel: ENSAMET
# btau: ENSBTAT
# cfam: ENSCAFT
# ecab: ENSECAT
# fcat: ENSFCAT
# hsap: ENST
# mfur: ENSMPUT
# elut: ELUT_
# pbra: PBRA_
sed -E -i'' 's/^>ENSAMET.*$/>amel/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ENSBTAT.*$/>btau/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ENSCAFT.*$/>cfam/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ENSECAT.*$/>ecab/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ENSFCAT.*$/>fcat/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ENST0.*$/>hsap/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna
sed -E -i'' 's/^>ENSMPUT.*$/>mfur/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ELUT_.*$/>elut/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
## add on PBRA:
sed -E -i'' 's/^>PBRA_.*$/>pbra/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna
# modified umar, lwed and oros to have SPP_ name in front of the protein name to make this more reliable
sed -E -i'' 's/^>UMAR_.*$/>umar/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna
sed -E -i'' 's/^>OROS_.*$/>oros/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna
sed -E -i'' 's/^>LWED_.*$/>lwed/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna
## add on mluc, pvam, ttru:
sed -E -i'' 's/^>ENSMLUT.*$/>mluc/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ENSPVAT.*$/>pvam/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 
sed -E -i'' 's/^>ENSTTRT.*$/>ttru/g' group_${j}/cluster_${cluster}_group_${j}_${TName}/Seqs_cluster_${cluster}_group_${j}_${TName}.fna 

done

sleep 5m
