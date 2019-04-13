## SIMPLE r script to select longest OR from results of step 3
# set minimum length to 40 AA s (120bp) (like in my pos selection analysis)
# 3 categories:
# close to ends of scaff
# close to NNNs
# neither

input1="./ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.containNNN.stranded.IDlist"
input2="./ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.within30ofEnds.stranded.IDlist"
input3="./ORF_longThan_3_bp_step_2_results.olfacUniqueBlastHits.notCloseToEnds.noNNN.stranded.IDlist"
selectLongestOR <- function(input){
  # read in names
  d <- read.table(input,header=F)
  d$name <- sapply(strsplit(as.character(d$V1),"1-"),"[",1)
  # this loses strand info but don't need it here?
  d$aa <- sapply(strsplit(as.character(d$V1),"1-"),"[",2)
  d$length <- sapply(strsplit(d$aa,"_"),"[",1)
  dagg <- aggregate(as.numeric(as.character(d$length)), by = list(d$name), max)
  # split on the aa size
  longest <- paste(dagg$Group.1,"1-",dagg$x,"_aa",sep="") # this needs a space if input has been through fasta_tool 
  longest <- gsub(">","",x = longest) # get rid of >
  return(longest)
}


longest1 <- selectLongestOR(input1)
write.table(longest1,paste(input1,"LONGEST",sep="."),row.names=F,quote=F,col.names = F)


longest2 <- selectLongestOR(input2)
write.table(longest2,paste(input2,"LONGEST",sep="."),row.names=F,quote=F,col.names = F)

longest3 <- selectLongestOR(input3)
write.table(longest3,paste(input3,"LONGEST",sep="."),row.names=F,quote=F,col.names = F)


# 459 (~ 80 are above 250 AAs, those must be the ones i visually discarded, 80 sounds about right for that. The rest are too short. (premature stop codons? check)
