setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut2/")
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
pullOut_LengthFilter <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/elut2/step_1_c_CleanedByEvalue_ABScript/elut2.pullOut.lengthFilter120.fullInfo.txt",sep="\t",header=T) # these don't include too short genes, but include the best long hit that you want to use. Some best long hits got swapped for shorterhits with higher evalues. need to deal with that 
pullOut_noLF <- read.table("step_1_c_CleanedByEvalue_ABScript_noLengthFilter/elut2.cleanedByEvalue.noLengthFilter.Pruned.all.pullOut.info.txt",header=T,sep="\t")
# add clade info if not there already

############################## Separate Clade Info ################
pullOut_noLF$clade <- sapply(strsplit(as.character(pullOut_noLF$Query), "_Clade"), "[", 2) # got this from: https://stackoverflow.com/questions/17215789/extract-a-substring-in-r-according-to-a-pattern
pullOut_noLF$ORspp <- sapply(strsplit(as.character(pullOut_noLF$Query), "_Clade"), "[", 1)
head(pullOut_noLF)
############# FIGURE OUT HITS THAT HAVE CHANGED ########
# NOTE THAT THIS INCLUDES THE 1-2 GENES THAT GOT SWITCHED DUE TO LEN 50 FILTER
# ARG.
missing <- pullOut_LengthFilter[!(pullOut_LengthFilter$bedname %in% pullOut_noLF$bedname),]
dim(missing) # 113 Hmm that's a lot. Only mising 1 functional.
pullOut_noLF$class <- "class" # standing in
dim(pullOut_LengthFilter)
dim(pullOut_noLF)
# for now: 
pullOut_all <- pullOut_noLF # this doesn't include the missing stuff! be ware!
############### GOING TO HAVE TO DEAL WITH THIS MORE FORMALLY I GUESS? #######

##################### FUNCTIONAL ORS ##############
functionalORs <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/elut2/step_9_FinalFunctionalList/elut2.classI.and.classII.final.ORs.passedInspection.txt",header=F) # this is the result of step 9 of functional analysis
dim(functionalORs) 
# the names of these ORs have the extra info of length and "aa" at the end:
# elut_ScbS9RH_67127_851235_852047_1-348_aa
# want it to be
#elut_ScbS9RH_67127_851235_852047_1 to match your blast hits
# so
functionalORs_names <- sapply(strsplit(as.character(functionalORs$V1),"-[0-9]+_aa",perl=T),"[",1)
############### MISSING FUNCTIONAL GENES / TO EXCLUDE ###########
### because 1 pgene should be excluded 
#from step 1c: 
# none!
exclude <- "elut2_KZ291759.1_41704788_41705006_-1"
# and one functional gene isn't in pullOut
#missing
# find if any are missing in this not length filtered pullout:
missingFunc <- functionalORs_names[!(functionalORs_names %in% pullOut_all$bedname)]
# get the missing functional gene from the Functional length filter Pullout:
#1  missing 
pullOut_allPlusMissingFunc <- rbind(pullOut_all[!(pullOut_all$bedname %in% exclude),],pullOut_LengthFilter[pullOut_LengthFilter$bedname==missingFunc,])
#####################################################

################# NON FUNCTIONAL ###############
################## a. TRUNCATED ###############
truncORs_NC <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut2/step_4_b_visualize/VisualizationNotes/elut2.Nmissing.CMissing.IDList.format.txt",header=F,strip.white=T) # this is list of truncated ORs. 
dim(truncORs_NC) #26
truncORs_N <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut2/step_4_b_visualize/VisualizationNotes/elut2.NMissing.CComplete.IDList.format.txt",header=F,strip.white=T)
dim(truncORs_N) # 32
truncORs_C <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/elut2/step_4_b_visualize/VisualizationNotes/elut2.Ncomplete.Cmissing.IDList.format.txt",header=F,strip.white=T)
dim(truncORs_C) # 12
# THIS IS NOT TRUNCATED: 
#none of these around:

allTrunc <- unique(rbind(truncORs_C,truncORs_N,truncORs_NC)) # may not quite equal total because one may have doubled up if the ORFs were equally long -- using unique to eliminate 
# it's okay if sizes are close but don't quite match
sum(duplicated(rbind(truncORs_C,truncORs_N,truncORs_NC))) #1 that's where the dup is -- it's okay.  
colnames(allTrunc) <- "bedname"
dim(allTrunc)
################# b. possible pseudogenes (gotta dig into more) ################
# these are genes that are in pullOut_all, but aren't in truncated or in functional category. 
pgenes <- pullOut_allPlusMissingFunc[!(pullOut_allPlusMissingFunc$bedname %in% allTrunc$bedname) & !(pullOut_allPlusMissingFunc$bedname %in% functionalORs_names),]
dim(pgenes) #  
dim(allTrunc) # 
length(unique(functionalORs_names)) # 571

### CHECK THIS CAREFULLY it was +1 because of hits that changed with and without length filter
### because 1 pgene should be excluded 
#from step 1c: 
#exclude <- "mfur_GL897276.1_595480_596211_-1"
# and one functional gene isn't in pullOut
#missing


############################## LABEL ALL GENES ###############
pullOut_allPlusMissingFunc$label <- NA
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$bedname %in% pgenes$bedname,]$label <- "pseudogene"
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$bedname %in% functionalORs_names,]$label <- "functional"
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$bedname %in% allTrunc$bedname,]$label <- "truncated"
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$bedname %in% pgenes$bedname & pullOut_allPlusMissingFunc$Hit_length < 200,]$label <- "fragment"
# check:
sum(pullOut_allPlusMissingFunc$label=="pseudogene") # 234 # why is this now 447, where is this extra 70 coming from?
sum(pullOut_allPlusMissingFunc$label=="functional") # 571
sum(pullOut_allPlusMissingFunc$label=="truncated") # 70 
sum(pullOut_allPlusMissingFunc$label=="fragment") # 353 



################# PLOT GENES BY CLADE ################
# # calculate % pseudogenized:
library(tidyr)
library(dplyr)
# fix awkward _CLASS_I
pullOut_allPlusMissingFunc$clade <- as.factor(pullOut_allPlusMissingFunc$clade)

levels(pullOut_allPlusMissingFunc$clade)[levels(pullOut_allPlusMissingFunc$clade)=="_CLASS_I"]<- "ClassI"
#pullOut_all[pullOut_all$clade=="_CLASS_I",]$clade <- "ClassI"
# https://stackoverflow.com/questions/40348854/geom-bar-percentage-for-each-group-in-x
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$clade=="ClassI",]$class <- "I"
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$clade!="ClassI",]$class <- "II"
## write out:
write.table(pullOut_allPlusMissingFunc,"step_5_labelAllgenes/pullOut.func.trunc.pgene.labelled.elut2.all.AddedMissing.Fragments.txt",row.names=F,quote = F,sep="\t")

tally1 <- pullOut_allPlusMissingFunc %>% group_by(clade,label) %>% tally %>% mutate(n/sum(n)) %>% mutate(sum(n))
# this is so cool! I need to spend more time with dplyr/tidyr
tally1
colnames(tally1) <- c("clade","label","n","perc","sumN")
head(tally1)

p2 <- ggplot(pullOut_allPlusMissingFunc,aes(x=clade,fill=label))+
  geom_bar(stat="count")+
  theme_bw()+
  scale_fill_manual(values=c(func,pg,tr))+
  ggtitle("Sea Otter (kenyoni) Olfactory Receptors")+
  theme(legend.title = element_blank(),legend.position="top",legend.background = element_rect("transparent"))+
  geom_text(data=tally1[tally1$label=="pseudogene",], aes(x=clade,y=sumN+10,label=n),size=4,nudge_x = 0,nudge_y=0.6,color=pg)+
  theme(axis.text.x = element_text(angle=45, hjust=1)) 
p2

ggsave("step_5_labelAllgenes/labeledGeneBreakdown.countPgeneLabel.pdf",p2,device="pdf",width=7,height=4)
