############ PACKAGES ###############
require(ape)
#install.packages('dendextend')
require(dendextend)
require(ggplot2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")
#biocLite("GenomeGraphs")
require(ggtree)
require(geiger)
require(GenomeGraphs)
############ WORKING DIRECTORY ############
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/pbra/")
###### READ IN TREE  ####################
#  this is from clustalX with distance correction and gaps removed
# and bootstraps as node labels
tree <- read.tree("step_6_alignWithOutgroups/step_6_result_mafftAligment.wOutgroups.20171003.2.phb")

############ SET UP GROUPS ###########
# outgroup
outgroup <- tree$tip.label[grep("pbra",tree$tip.label,invert=T)]
outgroup <- outgroup[outgroup!="Human_OR2J3"]
outgroup <- outgroup[grep("Clade",outgroup,invert=T)]
pbra <- tree$tip.label[grep("pbra",tree$tip.label,invert=F)]
# human
human <-  tree$tip.label[grep("Human",tree$tip.label,invert=F)]
# other spp:
representatives <- tree$tip.label[grep("Clade",tree$tip.label)]
classI <- representatives[grep("CLASS_I",representatives)]
classII <- representatives[grep("CLASS_I",representatives,invert=T)]
####### ROOT TREE USING OUTGROUP ###########
tree <- root(tree,outgroup = outgroup,resolve.root = T)

####### MAKE GROUPS IN THE TREE ########
# want to make each clade a group somehow?
representatives <- tree$tip.label[grep("Clade",tree$tip.label)]
clades <- sapply(strsplit(as.character(representatives), "_Clade"), "[", 2)

groups <- list(g1=outgroup,g2=pbra,g4=classI,g5=c(classII,human)) # create a list of your different groups! nice!

# to color a group
tree <- groupOTU(tree,groups)

######################## PLOT #########################
p <- ggtree(tree,layout="circular",aes(color=group))
p <- p + ggtitle("Giant Otter \n Class I and II Olfactory Receptors, with outgroup")+
  geom_tiplab(size=1,aes(angle=angle))+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3,size=0.5,color="black")
p

############# COLLAPSE CLADES  ############
# these numbers will change!!! 
cp <- p %>% collapse(node=577) %>% collapse(node=589)
cp <- cp + geom_point2(aes(subset=(node == 577)), size=5, shape=23, fill="steelblue")+
  geom_point2(aes(subset=(node == 589)), size=5, shape=23, fill="steelblue")
cp
################## OUTLIERs ###############
# seems to be one in class i 
viewClade(ggtree(tree)+geom_tiplab(),node=1107) #
outlier <- "pbra_flattened_line_22570_62447_63259_-1-348_aa"
tree.drop <- drop.tip(tree,outlier)
tree.drop <- groupOTU(tree.drop,groups)
p <- ggtree(tree.drop,layout="circular",aes(color=group)) + ggtitle("Giant Otter \n Class I and II Olfactory Receptors, with outgroup")+
geom_tiplab(size=1,aes(angle=angle))+
geom_text2(aes(subset=!isTip, label=node), hjust=-.3,size=0.5,color="black")
p
# class I: node 1794
# class II: 2 nodes (1906, 976)
############### Figure out Class I and Class II nodes ###########
# class I:  node
classINode <- 1794
viewClade(ggtree(tree.drop)+geom_tiplab(size=3), node=classINode)
classITips <- tips(tree.drop,node=classINode)
classITips_sppOnly <- classITips[grep("pbra",classITips)] #  
length(classITips_sppOnly) # 101 (77 in elut)
# class II;  node 1068 and 1079 (clade B)
classIINode.1 <- 1906
classIINode.2 <- 976
viewClade(ggtree(tree.drop)+geom_tiplab(), node=classIINode.1)
classIITips.1 <- tips(tree.drop,node=classIINode.1)
classIITips.2 <- tips(tree.drop,node=classIINode.2)
classIITips <- c(classIITips.1,classIITips.2)
classIITips_sppOnly <- classIITips[grep("pbra",classIITips)] # 502 
length(classIITips_sppOnly)
# these contain the non-pbra genes as well.

length(classITips)
length(classIITips)

# # write out:
 write.table(classITips,"step_8_alignClassesSeparately/classI.tips.includingReps.noOutlier.txt",row.names=F,col.names=F,quote=F)
 write.table(classITips_sppOnly,"step_8_alignClassesSeparately/classI.tips.sppOnly.noOutlier.txt",row.names=F,col.names=F,quote=F)
 write.table(classIITips,"step_8_alignClassesSeparately/classII.tips.includingReps.noOutlier.txt",row.names=F,col.names=F,quote=F)
 write.table(classIITips_sppOnly,"step_8_alignClassesSeparately/classII.tips.sppOnly.noOutlier.txt",row.names=F,col.names=F,quote=F)
############ ZOOM makes nice figure ##############
gzoom(tree,classITips)
###### Plot positions 
classIlocations <- strsplit(classITips_sppOnly,"_")

classIscaffs <- sapply(classIlocations, "[", 3)
classIstarts <- sapply(classIlocations, "[", 4)
classIscaff_start <- as.data.frame(cbind(classIscaffs,classIstarts))


classIweirdclade <- strsplit(tips(tree,1290),"_")
classIweirdclade_starts <- sapply(classIweirdclade, "[", 4)
classIweirdclade_scaffs <- sapply(classIweirdclade, "[", 3)
classIweirdclade_scaff_start <- as.data.frame(cbind(classIweirdclade_scaffs,classIweirdclade_starts))


ggplot(classIscaff_start,aes(y=as.character(classIscaffs),x=as.numeric(as.character(classIstarts))))+
  geom_point()+
  geom_point(data=classIweirdclade_scaff_start,color="pink",aes(y=as.character(classIweirdclade_scaffs),x=as.numeric(as.character(classIweirdclade_starts))))

classIIlocations <- strsplit(classIITips_sppOnly,"_")
classIIlocations <- as.table(classIIlocations)
classIIscaffs <- sapply(classIIlocations, "[", 3)
classIIstarts <- sapply(classIIlocations, "[", 4)
classIIscaff_start <- as.data.frame(cbind(classIIscaffs,classIIstarts))

ggplot(classIIscaff_start,aes(y=as.character(classIIscaffs),x=as.numeric(as.character(classIIstarts))))+
  geom_point()
