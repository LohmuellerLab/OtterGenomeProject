### Idea:
# tile the genome into even windows
# compute ABBA BABA counts into those windows
# then jacknife (calculate D stat dropping each one of those windows at a time)
data.dir="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/resultsCombining_all3Otters_inclNSO/ABBA-BABA/"
# set bin size you want (want it to be longer than recombination block)
binsize=500000 # 500kb
# to get tiles:
#https://bioinformatics.stackexchange.com/questions/2924/dividing-genome-into-non-overlapping-windows-using-r
require(GenomicRanges)
require(bootstrap)
# to install genomicRanges:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenomicRanges", version = "3.8")
# if you need to do from scratch

################## sum abba-baba within ranges #########################
######## TEMPORARY!!!! 
abbaBaba <- read.table(paste(data.dir,"results/ABBA-BABA.allScaffs.txt.gz",sep=""),header=T,sep="\t",stringsAsFactors = F)
head(abbaBaba)
# incorporate into GRanges:
# make Granges object of AbbaBaba:
# another useful function:makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
abbaBaba_GRanges <- GRanges(seqnames=as.character(abbaBaba$scaffold),ranges=as.character(abbaBaba$pos),abba=as.numeric(abbaBaba$abba),baba=as.numeric(abbaBaba$baba))
abbaBaba_GRanges

####################### set up bins ################
# get the names and lengths of all chromosomes (scaffolds)
# I had already gotten chr sizes for work in the past from Mustela genome:
mustelaChrSizes <- (read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/SlidingWindowheterozygosity/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.224.ChrLengths.txt",sep=",",stringsAsFactors = F,strip.white = T))
mustelaChrSizes <- data.frame(t(mustelaChrSizes))
colnames(mustelaChrSizes) <- "record"
mustelaChrSizes$scaff <- lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",1)
mustelaChrSizes$size <- lapply(strsplit(as.character(mustelaChrSizes$record),":"),"[",2)
mustelaChrSizes_AB <-  mustelaChrSizes[mustelaChrSizes$scaff %in% seqnames(abbaBaba_GRanges),] # select those for which we have abba baba information
# only get sizes for which you have abba baba sites; rest are uniformative and will break things
sizes <- as.numeric(mustelaChrSizes_AB$size)
names(sizes) <- mustelaChrSizes_AB$scaff
head(sizes)
# 
#     Chr1     Chr2     Chr3     Chr4     Chr5 
# 18585056 19698289 23459830 26975502 30427671 
#  note you can also make sliding windows this way! 

### 500kb bins for jackknifing
bins   <- tileGenome(sizes, tilewidth=binsize, cut.last.tile.in.chrom=T)
bins
# add a bin number to each bin
bins$binNum <- seq(1,length(bins))
length(unique(bins$binNum))

### make RLEs :

abbaRLE <- mcolAsRleList(abbaBaba_GRanges, varname = "abba")
babaRLE <- mcolAsRleList(abbaBaba_GRanges, varname = "baba")
# function to count items per bin from: http://crazyhottommy.blogspot.com/2016/02/compute-averagessums-on-granges-or.html
# added na.rm=T
binnedSum <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewSums(views,na.rm = T)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

# bins need to just apply to abbaRLE and make it a dataframe:
abba_sum_binned = as.data.frame(binnedSum(bins,abbaRLE,"abba_sum"))
baba_sum_binned = as.data.frame(binnedSum(bins,babaRLE,"baba_sum"))

# omg this worked!!!

# merge together and get ready for jackknifing:
abba_baba_bin_summed <- merge(abba_sum_binned[,c("binNum","abba_sum")],baba_sum_binned[,c("binNum","baba_sum")],by="binNum")

head(abba_baba_bin_summed)

# okay I can now jackknife on this!! 
#https://www.rdocumentation.org/packages/bootstrap/versions/2015.2/topics/jackknife
# define function -- this is defined in a really WEIRD way but it works
# x is NOT a row number, but is the set of all rows that you want will be INCLUDED in the set
# so sum(dataframe[x,]$abba_sum) will be the whole set in minus one row
# then the data you pass to jackknife is just a vector of numbers that are the length of your dataframe
# and this tells it to go through that many jackknifes. where the number x represents the WHOLE SET minus one 
# this is totally counterintutive and weird. But I tested with a dummy dataset and it works perfectly if you do it this way. Confirmed with this blog post as well https://andyphilips.github.io/blog/2017/01/09/bootstrap-jackknife.html 
# x is the 
Dstat_binNum=function(x,dataframe){
  # abba_baba_bin_summed is in format binNum abba_sum baba_sum where abba and baba sums are sums of that bin
  D=(sum(dataframe[x,]$abba_sum)-sum(dataframe[x,]$baba_sum))/(sum(dataframe[x,]$abba_sum)+sum(dataframe[x,]$baba_sum))
  return(D)
}

jack <- jackknife(1:length(abba_baba_bin_summed$binNum),Dstat_binNum,abba_baba_bin_summed) # this worked!
jack$jack.se # standard error
# from http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/
D_noJack <- (sum(abba_baba_bin_summed$abba_sum) - sum(abba_baba_bin_summed$baba_sum))/(sum(abba_baba_bin_summed$abba_sum) + sum(abba_baba_bin_summed$baba_sum)) # actual estimate of D from all the data
D_noJack
D_err <- jack$jack.se # standard error 
D_Z <- D_noJack / D_err # z-score (x - mean / sterror); here mean is 0 because we are testing if ABBA-BABA is sig diff from zero. so it just becomes D / DstdErr
D_p <- 2*pnorm(-abs(D_Z)) # pvalue
D_p
sink(paste(data.dir,"results/ABBA-BABA.Jackknife.txt",sep=""))
cat("D-Stat from data: ",D_noJack,"\n")
cat("Jackknife std error: ",D_err,"\n")
cat("D Z score (D/D_stdError): ",D_Z,"\n") 
cat("D p value: ",D_p,"\n")
sink()
