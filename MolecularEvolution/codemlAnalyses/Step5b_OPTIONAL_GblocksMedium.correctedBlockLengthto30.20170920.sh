#! /bin/bash
#$ -l h_rt=5:00:00,h_data=10G
#$ -cwd
#$ -N gblocksTest
#$ -m abe
#$ -M ab08028
#$ -o gblocksReports
#$ -e gblocksReports
#$ -t 1-77
gblocks=/u/home/a/ab08028/klohmueldata/annabel_data/bin/Gblocks

# this will convert output to phylip format (from: https://indra.mullins.microbiol.washington.edu/cgi-bin/perlscript/info.cgi?ID=Fasta2Phylip.pl&path=perlscript-scripts)

convert=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/scripts/fasta2phylip.pl

res=0.93 # residue masking filter I used 
row1=11 # row where hsap ID is
row2=8 # alt name if no hsap ID, use can fam
clusters=one2one_orthologClusters_wideFormat_15spp.includesNAs.20170606.txt # clusters file 
date=20170608 # date you ran guidance
dirDate=20170606
j=`printf "%02d" $((10#$SGE_TASK_ID - 10#1))`
start=$((10#$j*200))
echo "start: " $start
end=$((10#$j*200+199))
echo "end: " $end
for (( i=$start; i<=$end; i++ ))
do 
echo $i
p=$(($i+1))
# get human T ID to be overall gene name
TName=`awk -v a=${row1} 'FNR==a {print $"'$p'"}' $clusters`
# but if TName is missing in human, sub in the dog:
if [[ $TName = "NA" ]]
then 
TName=`awk -v a=${row2} 'FNR==a {print $"'$p'"}' $clusters`
fi

cluster=`printf "%05d" $i`
wd=/u/home/a/ab08028/klohmueldata/annabel_data/comparativeGenomicsDec2016/sequenceAlignment_${dirDate}
cd $wd
header=group_${j}/cluster_${cluster}_group_${j}_${TName}
guidanceRun=guidance_${date}
fastaFile=$header/$guidanceRun/MSA.PRANK.aln.With_Names.ResMask.${res}

if [ -e $fastaFile ]
then
echo found $fastaFile
# get spp count:
sppCount=`grep -c "^>" $fastaFile`
b2_value=$((((sppCount+1)/2)+1)) 
echo b2-value: $b2_value
# b2_value: this rounds up halves (want half plus 1, but need to make an integer; so test=17 would give (18/2)+1 = 10, which is rounded up from 9.5)
############### GBLOCKS parameters
# b4: min block size (in nucleotides)
# t : codons
# b5: half of gaps 
# b3: max contiugous conserved (default = 8; relaxed =10)
# b2: min for flank (half plus 1)
sed 's/n/N/g' $fastaFile > ${fastaFile}.n2N.forGblocks.fasta
# 20170920  changing blcoks block size to 30! (nucleotides not codons) 
# adding e: output suffix gbmed (must be 5 chars)
$gblocks ${fastaFile}.n2N.forGblocks.fasta -t=c -b5=h -b4=30 -b3=10 -b2=${b2_value} -e=gbmed
# note: gblocks doesn't like paths that start with /u/home/.... but can deal with being at the higher level
# as long as it's relative to where you are 
sed -E 's/([ATGCN\-])( )([NATGC\-])/\1\3/g' ${fastaFile}.n2N.forGblocks.fastagbmed > ${fastaFile}gbmed.NOSPACES

# this makes it phylip format:
$convert ${fastaFile}gbmed.NOSPACES > ${fastaFile}gbmed.phylip
fi
done

sleep 10m
