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
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/elut1/elut1.duplicates/")
###### READ IN TREE  ####################
#  this is from clustalX with distance correction and gaps removed
# and bootstraps as node labels
tree <- read.tree("step_6_MakeNJTree/step_6_result_mafftAligment.wOutgroups.20171103.phb")

############ SET UP GROUPS ###########
# outgroup
# pbra --> pbra
outgroup <- tree$tip.label[grep("elut",tree$tip.label,invert=T)]
outgroup <- outgroup[outgroup!="Human_OR2J3"]
outgroup <- outgroup[grep("Clade",outgroup,invert=T)]
# sea otter (eventually do this for pbra and mfur?)
elut<- tree$tip.label[grep("elut",tree$tip.label,invert=F)]
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
  ggtitle("Sea Otter Dup \nClass I and II Olfactory Receptors, with outgroup")+
  geom_tiplab(aes(label=label2),size=1.5,alpha=0.8,color="black")+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3,size=0.5,color="black")
# to get bootstrap labels, label=label
# to get node labels, label=node
# shorten labels? rememver geom_label_repel for later
p
# sweet all look good. 
#p %<+% d # weird operator from ggtree, attaches the info on based on first column being the same 
#https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html#the-operator

############### Figure out Class I and Class II nodes ###########
# can't resolve trichotomy? more bootstraps? different program?
# class I:  
classInode = 749
viewClade(ggtree(tree)+geom_tiplab(size=3), node=classInode)
classITips <- tips(tree,node=classInode)
classITips_sppOnly <- classITips[grep("elut",classITips)] 
length(classITips_sppOnly) # +21 dup (77 in elut) --> 141 in mfur (incomplete genome or pseudogenes?)
# class II;  
classIInode  = 738
viewClade(ggtree(tree)+geom_tiplab(), node=classIInode)
classIITips <- tips(tree,node=classIInode)
classIITips_sppOnly <- classIITips[grep("elut",classIITips)] # (elut: 398), mfur: 682
length(classIITips_sppOnly) # 42 new from dup. 
# these contain the non-pbra genes as well.
# total of 63 new seqs
length(classITips)
length(classIITips)
############ REMOVE OUTLIER(s) ############################
# write out:
write.table(classITips,"step_8_alignClassesSeparately/classI.tips.includingReps.noOutlier.txt",row.names=F,col.names=F,quote=F)
write.table(classITips_sppOnly,"step_8_alignClassesSeparately/classI.tips.sppOnly.noOutlier.txt",row.names=F,col.names=F,quote=F)
write.table(classIITips,"step_8_alignClassesSeparately/classII.tips.includingReps.noOutlier.txt",row.names=F,col.names=F,quote=F)
write.table(classIITips_sppOnly,"step_8_alignClassesSeparately/classII.tips.sppOnly.noOutlier.txt",row.names=F,col.names=F,quote=F)

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
