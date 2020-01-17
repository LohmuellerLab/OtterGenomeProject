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
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults")
###### READ IN TREE  ####################
#  this is from clustalX with distance correction and gaps removed
# and bootstraps as node labels
tree <- read.tree("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/step_6_result_mafftAligment.wOutgroups.wRepresentatives.noOutlier.20170929.phb")

############ SET UP GROUPS ###########
# outgroup
outgroup <- tree$tip.label[grep("elut",tree$tip.label,invert=T)]
outgroup <- outgroup[outgroup!="Human_OR2J3"]
outgroup <- outgroup[grep("Clade",outgroup,invert=T)]
# sea otter (eventually do this for pbra and mfur?)
elut <- tree$tip.label[grep("elut",tree$tip.label,invert=F)]
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

groups <- list(g1=outgroup,g2=elut,g4=classI,g5=c(classII,human)) # create a list of your different groups! nice!

# to color a group
tree <- groupOTU(tree,groups)

################# USEFUL THINGS: MRCA, viewClade, etc. ########
# to get a node #:
MRCA(tree, tip=outgroup) # 579 
MRCA(tree,tip=elut) # 489
 # node 579 highligh
viewClade(ggtree(tree)+geom_tiplab(), node=579) # aha! this zooms in on a node
viewClade(ggtree(tree)+geom_tiplab(), node=589) # aha! this zooms in on a node

viewClade(ggtree(tree)+geom_tiplab(), node=489) # aha! this zooms in on a node I think this is Class I!
viewClade(ggtree(tree)+geom_tiplab(), node=576)

######################## PLOT #########################
p <- ggtree(tree,layout="circular",aes(color=group))
p <- p + ggtitle("Class I and II Olfactory Receptors, with outgroup")+
  geom_tiplab(size=1,aes(angle=angle))+
  geom_text2(aes(subset=!isTip, label=label), hjust=-.3,size=0.5,color="black")
p

############# COLLAPSE CLADES  ############
# these numbers will change!!! 
cp <- p %>% collapse(node=577) %>% collapse(node=589)
cp <- cp + geom_point2(aes(subset=(node == 577)), size=5, shape=23, fill="steelblue")+
  geom_point2(aes(subset=(node == 589)), size=5, shape=23, fill="steelblue")
cp
############### BOOTSTRAP VALUES #####################
############### Figure out Class I and Class II nodes ###########
# class I: 1227 node
viewClade(ggtree(tree)+geom_tiplab(size=3), node=1227)
classITips <- tips(tree,node=1227)
classITips_sppOnly <- classITips[grep("elut",classITips)] # 77 
# class II; 1225 node
viewClade(ggtree(tree)+geom_tiplab(), node=1225)
classIITips <- tips(tree,node=1225)
classIITips_sppOnly <- classIITips[grep("elut",classIITips)] # 398 
length(classIITips_sppOnly)
# these contain the non-elut genes as well.

length(classITips)
length(classIITips)

# write out:
write.table(classITips,"listsOfORs/elut/classI.tips.includingReps.txt",row.names=F,col.names=F,quote=F)
write.table(classITips_sppOnly,"listsOfORs/elut/classI.tips.sppOnly.txt",row.names=F,col.names=F,quote=F)
write.table(classIITips,"listsOfORs/elut/classII.tips.includingReps.txt",row.names=F,col.names=F,quote=F)
write.table(classIITips_sppOnly,"listsOfORs/elut/classII.tips.sppOnly.txt",row.names=F,col.names=F,quote=F)
############ ZOOM makes nice figure ##############
gzoom(tree,classITips)
###### Plot positions  #########
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
