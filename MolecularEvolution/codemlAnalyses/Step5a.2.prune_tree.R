args <- commandArgs(trailingOnly = TRUE)
# usage: Step5a.2.prune_tree.R <mainTree> <header>
library(ape)
tree=read.tree(args[1])
taxa = readLines(paste(args[2],"/taxaList.txt",sep=""))
pruned.tree <- drop.tip(tree, setdiff(tree$tip.label,taxa))
write.tree(pruned.tree,file=paste(args[2],"/prunedTree.noFG.txt",sep=""))
