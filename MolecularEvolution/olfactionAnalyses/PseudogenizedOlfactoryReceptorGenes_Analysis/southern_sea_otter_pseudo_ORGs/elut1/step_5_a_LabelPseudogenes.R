setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut1/")
require(ggplot2)
require(RColorBrewer)
# colors:
#RColorBrewer::display.brewer.all()
RColorBrewer::display.brewer.pal(name="Set1",n=3)
cols <- brewer.pal(name="Set1",n=3)
tr <- cols[2]
func <- cols[3]
pg <- cols[1]
########### ALL HITS, not length filtered ############# 
pullOut_all <- read.table("step_1_c_CleanedByEvalue_ABScript_noLengthFilter/pullOut.cleanedByEvalue.noLengthFilter.Pruned.allHits.FuncAndNonFunc.txt",header=T,sep="\t")
##################### FUNCTIONAL ORS ##############
functionalORs <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/elut1/step_9_FinalFunctionalList/classI.and.classII.final.ORs.passedInspection.txt",header=F) # this is the result of step 9 of functional analysis
dim(functionalORs) 
# the names of these ORs have the extra info of length and "aa" at the end:
# elut_ScbS9RH_67127_851235_852047_1-348_aa
# want it to be
#elut_ScbS9RH_67127_851235_852047_1 to match your blast hits
# so
functionalORs_names <- sapply(strsplit(as.character(functionalORs$V1),"-[0-9]+_aa",perl=T),"[",1)


#####################################################
############# FIGURE OUT HITS THAT HAVE CHANGED ########
functionalORs_names[!(functionalORs_names %in% pullOut_all$bedname)]
# cool, they're all there, don't need to do missingness thing.
################# NON FUNCTIONAL ###############
################## a. TRUNCATED ###############
truncORs_NC <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut1/step_4_b_visualize/VisualizationNotes/elut1.Nmissing.Cmissing.IDList.format.txt",header=F,strip.white=T) # this is list of truncated ORs. 
truncORs_N <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut1/step_4_b_visualize/VisualizationNotes/elut1.Nmissing.Ccomplete.IDList.format.txt",header=F,strip.white=T)
truncORs_C <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut1/step_4_b_visualize/VisualizationNotes/elut1.Ncomplete.Cmissing.IDList.format.txt",header=F,strip.white=T)
notTrunc <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut1/step_4_b_visualize/VisualizationNotes/Ncomplete.Ccomplete.Pgene.IDList.format.txt",header=F,strip.white = T)
allTrunc <- unique(rbind(truncORs_C,truncORs_N,truncORs_NC)) # one may have doubled up
colnames(allTrunc) <- "bedname"
################# b. possible pseudogenes (gotta dig into more) ################
# these are genes that are in pullOut_all, but aren't in truncated or in functional category. 
pgenes <- pullOut_all[!(pullOut_all$bedname %in% allTrunc$bedname) & !(pullOut_all$bedname %in% functionalORs_names),]
dim(pgenes)
dim(allTrunc)
length(functionalORs_names)
# 461+185+468 = 1114 ! this is correct.

############################## LABEL ALL GENES ###############
pullOut_all$label <- NA
pullOut_all[pullOut_all$bedname %in% pgenes$bedname,]$label <- "pseudogene"
pullOut_all[pullOut_all$bedname %in% functionalORs_names,]$label <- "functional"
pullOut_all[pullOut_all$bedname %in% allTrunc$bedname,]$label <- "truncated"
pullOut_all[pullOut_all$bedname %in% pgenes$bedname & pullOut_all$Hit_length < 200,]$label <- "fragment"
# check:
sum(pullOut_all$label=="pseudogene") # 183
sum(pullOut_all$label=="functional") #  468
sum(pullOut_all$label=="truncated") # 185
sum(pullOut_all$label=="fragment") # 278
# cool. 
write.table(pullOut_all,"step_5_labelAllgenes/pullOut.func.trunc.pgene.labelled.elut1.all.fragments.txt",row.names=F,quote = F,sep="\t")
################# PLOT GENES BY CLADE ################
# fix awkward _CLASS_I

levels(pullOut_all$clade)[levels(pullOut_all$clade)=="_CLASS_I"]  <- "ClassI"

# # calculate % pseudogenized:
library(tidyr)
library(dplyr)
# https://stackoverflow.com/questions/40348854/geom-bar-percentage-for-each-group-in-x
tally1 <- pullOut_all %>% group_by(clade,label) %>% tally %>% mutate(n/sum(n)) %>% mutate(sum(n))
# this is so cool! I need to spend more time with dplyr/tidyr
tally1
colnames(tally1) <- c("clade","label","n","perc","sumN")
head(tally1)
write.table(tally1,"step_5_labelAllgenes/tally.byCladeCategory.txt",row.names = F,quote=F,sep="\t")
# fix awkward _CLASS_I
# plot:
p <- ggplot(pullOut_all,aes(x=clade,fill=label))+
  geom_bar(stat="count")+
  theme_bw()+
  scale_fill_manual(values=c(func,pg,tr))+
  ggtitle("Sea Otter (nereis) Olfactory Receptors")+
  theme(legend.title = element_blank(),legend.position="top",legend.background = element_rect("transparent"))+
  geom_text(data=tally1[tally1$label=="pseudogene",], aes(x=clade,y=sumN+10,label=n),size=4,nudge_x = 0,nudge_y=0.2,color=pg)+
 theme(axis.text.x = element_text(angle=45, hjust=1)) 
p
ggsave("step_5_labelAllgenes/labeledGeneBreakdown.countPgeneLabel.pdf",p,device="pdf",width=7,height=4)

############ Class I and II ########################
pullOut_all$class <- NA
pullOut_all[pullOut_all$clade=="ClassI",]$class <- "I"
pullOut_all[pullOut_all$clade!="ClassI",]$class <- "II"
# write out with class info:
write.table(pullOut_all,"step_5_labelAllgenes/pullOut.func.trunc.pgene.labelled.elut1.all.fragments.txt",row.names=F,quote = F,sep="\t")
tally2 <- pullOut_all %>% group_by(class,label) %>% tally %>% mutate(n/sum(n)) %>% mutate(sum(n))
# this is so cool! I need to spend more time with dplyr/tidyr
tally2
colnames(tally2) <- c("class","label","n","perc","sumN")
head(tally2)
############# Plot Class I and II together #########
p <- ggplot(pullOut_all,aes(x=class,fill=label))+
  geom_bar(stat="count")+
  theme_bw()+
  scale_fill_manual(values=c(func,pg,tr))+
  ggtitle("Sea Otter (nereis) Olfactory Receptors")+
  theme(legend.title = element_blank(),legend.position="top",legend.background = element_rect("transparent"))+
  theme(axis.text.x = element_text(hjust=1)) 
p
################# Sandbox: Explore results ############
# I want to know a couple things:
# length distribution of pseudogene hits (missing some truncated ones?)
# how many of pgenes are class I vs Class I
head(pullOut_all)
ggplot(pullOut_all[pullOut_all$label=="pseudogene",],aes(x=Hit_length,fill=class))+
  geom_density(alpha=0.5)

ggplot(pullOut_all[pullOut_all$label=="truncated",],aes(x=Hit_length,fill=class))+
  geom_density(alpha=0.5)
