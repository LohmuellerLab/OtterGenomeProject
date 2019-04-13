
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

DP_data = read.table(args[1])
chrom = args[2]
outfile = args[3]

DP_sum = sum(as.numeric(DP_data[,1]))
DP_length = length(DP_data[,1])

toprint = data.frame(chrom, DP_sum, DP_length)

write.table(toprint, outfile, quote = FALSE, col.names = FALSE, row.names=FALSE)
