module load R
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/mapReads_General/northern_sea_otter_mapping
# script dir okay
indir=/u/flashscratch/a/ab08028/otters/vcfs/vcf_20190215_nso_lib283only/stats
for i in {1..50}
do
Rscript $scriptdir/compute_DP_sum.R $indir/DP.dist.scaff.${i}.out group_${i} DP.dist.scaff.${i}.sum
done 


# concatenate sums:
> $indir/50Scaffs_DP_SumsPerScaff.txt  
for i in {1..50}
do
cat DP.dist.scaff.${i}.sum >> $indir/50Scaffs_DP_SumsPerScaff.txt  
done
# then gzip the outs for ease of storage
gzip $indir/DP.dist.scaff.*.out