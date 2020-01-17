############ PACKAGES ###############
require(ape)
#install.packages('dendextend')
require(dendextend)
require(ggplot2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")
#biocLite("GenomeGraphs")
#biocLite("ggrepel")
require(ggtree)
require(geiger)
require(GenomeGraphs)
require(RColorBrewer)
require(ggrepel)
############ WORKING DIRECTORY ############
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/elut2/")
###### READ IN TREE  ####################
#  this is from clustalX with distance correction and gaps removed
# and bootstraps as node labels
tree <- read.tree("step_6_alignWithOutgroups/step_6_result_mafftAligment.wOutgroups.20171003.phb")

############ SET UP GROUPS ###########
# outgroup

outgroup <- tree$tip.label[grep("elut2",tree$tip.label,invert=T)]
outgroup <- outgroup[outgroup!="Human_OR2J3"]
outgroup <- outgroup[grep("Clade",outgroup,invert=T)]
elut2 <- tree$tip.label[grep("elut2",tree$tip.label,invert=F)]
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


############# SHORTEN CLADE LABELS ##############
origlabel <- tree$tip.label
d <- data.frame(origlabel=origlabel,label2=NA)
# name the clades:
d[grep("Clade",d$origlabel),]$label2 <-  sapply(strsplit(as.character(d[grep("Clade",d$origlabel),]$origlabel),"_Clade"),"[",2)
# label class I:
d[grep("CLASS_I",d$origlabel),]$label2 <-  "Class I"
######################## PLOT #########################
p <- ggtree(tree,layout="circular",aes(color=group)) 
p <- p %<+% d # attaches new labels
p <- p + 
  ggtitle("Sea Otter 2 (kenyoni) \nClass I and II Olfactory Receptors, with outgroup")+
  geom_tiplab(aes(label=label2),size=1.5,alpha=0.8,color="black")+
  geom_text2(aes(subset=!isTip, label=label), hjust=-.3,size=0.5,color="black")
# to get bootstrap labels, label=label
# to get node labels, label=node
p
#p %<+% d # weird operator from ggtree, attaches the info on based on first column being the same 
#https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html#the-operator
########### OUTLIER ############
# look with node #s:
p <- ggtree(tree,layout="circular",aes(color=group)) 
p <- p %<+% d # attaches new labels
p <- p + 
  ggtitle("Ferret \nClass I and II Olfactory Receptors, with outgroup")+
  geom_tiplab(aes(label=label2),size=1.5,alpha=0.8,color="black")+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3,size=0.5,color="black")
# to get bootstrap labels, label=label
# to get node labels, label=node
p
outlier <- "elut2_GL897417.1_82828_83739_1-323_aa" # visual inspection. Remove this sequence
outlierNode=1679
viewClade(ggtree(tree)+geom_tiplab(size=.5),node=outlierNode)
############# COLLAPSE CLADES  ############
cp <- p %>% collapse(node=577) %>% collapse(node=589)
cp <- cp + geom_point2(aes(subset=(node == 577)), size=5, shape=23, fill="steelblue")+
  geom_point2(aes(subset=(node == 589)), size=5, shape=23, fill="steelblue")
cp

############### Figure out Class I and Class II nodes ###########
# class I:
classInode = 1689
viewClade(ggtree(tree)+geom_tiplab(size=3), node=classInode)
classITips <- tips(tree,node=classInode)
classITips_sppOnly <- classITips[grep("elut2",classITips)] 
length(classITips_sppOnly) # (77 in elut) --> 141 in elut2 (incomplete genome or pseudogenes?)
# class II;  
classIInode  = 1677
viewClade(ggtree(tree)+geom_tiplab(), node=classIInode)
classIITips <- tips(tree,node=classIInode)
classIITips_sppOnly <- classIITips[grep("elut2",classIITips)] # (elut: 398), elut2: 682
length(classIITips_sppOnly)
# these contain the non-pbra genes as well.

length(classITips)
length(classIITips)
############ REMOVE OUTLIER(s) ############################
classITips <- classITips[!(classITips %in% outlier)]
classIITips <- classIITips[!(classIITips %in% outlier)]

# write out:
 write.table(classITips,"step_8_alignClassesSeparately/classI.tips.includingReps.noOutlier.txt",row.names=F,col.names=F,quote=F)
 write.table(classITips_sppOnly,"step_8_alignClassesSeparately/classI.tips.sppOnly.noOutlier.txt",row.names=F,col.names=F,quote=F)
 write.table(classIITips,"step_8_alignClassesSeparately/classII.tips.includingReps.noOutlier.txt",row.names=F,col.names=F,quote=F)
 write.table(classIITips_sppOnly,"step_8_alignClassesSeparately/classII.tips.sppOnly.noOutlier.txt",row.names=F,col.names=F,quote=F)
############ ZOOM makes nice figure ##############
gzoom(tree,classITips)

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
