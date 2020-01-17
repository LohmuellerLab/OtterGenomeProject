# check that all 100%
# check that name matches closely (w/in some amount?)
# make final bed file . need to treat + and - strands differently
# usage: Rscript <script.R> spp
args <- commandArgs()
spp=args[1]
spp="elut2"
setwd(paste("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/",spp,sep=""))
input <- read.table("step_10_getFinalGenomeCoordsAndSeqs/bigger.test.blastOutput.txt",header=F) # this is the result of step 10a
# : tblastn -query $query -db $db -num_threads 20 -outfmt 6 -max_target_seqs 1 > step_10_getFinalGenomeCoordsAndSeqs/step_10_finalFunctionalPGenes_blastOutput.txt
# 
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
colnames(input) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
head(input)
# check that all are 100 percent
input[input$pident != 100,]

# check that name matches
input$origScaff <- lapply(strsplit(as.character(input$qseqid),"_"),"[",2)
input[input$origScaff != input$sseqid,]

# get strand myself:
input$strand <- NA
input[as.numeric(input$send) > as.numeric(input$sstart),]$strand <- "+"
input[as.numeric(input$send) < as.numeric(input$sstart),]$strand <- "-"
# (check it)
View(input)
# make bed:
input$bedStart <- NA
input$bedEnd <- NA
input[input$strand=="+",]$bedStart <- input[input$strand=="+",]$sstart - 1
input[input$strand=="-",]$bedStart <- input[input$strand=="-",]$send - 1
input[input$strand=="+",]$bedEnd <- input[input$strand=="+",]$send
input[input$strand=="-",]$bedEnd <- input[input$strand=="-",]$sstart

bed <- input[,c("sseqid","bedStart","bedEnd","qseqid","evalue","strand")]

write.table(bed,paste("step_10_getFinalGenomeCoordsAndSeqs/",spp,".FinalFunctionalORs.coords.TEST.bed",sep=""),row.names = F,col.names = F,quote=F,sep="\t")


# write out coordinates in step10
