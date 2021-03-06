setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/mfur/")
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
pullOut_LengthFilter <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/mfur/step_1_c_CleanedByEvalue_ABScript/pullOut.mfur.fullInfo.Lengthrestriction250.txt",sep="\t",header=T) # these don't include too short genes, but include the best long hit that you want to use. Some best long hits got swapped for shorterhits with higher evalues. need to deal with that 
pullOut_noLF <- read.table("step_1_c_CleanedByEvalue_ABScript_noLengthFilter/mfur.cleanedByEvalue.noLengthFilter.Pruned.all.pullOut.info.txt",header=T,sep="\t")
# add clade info if not there already

############################## Separate Clade Info ################
pullOut_noLF$clade <- sapply(strsplit(as.character(pullOut_noLF$Query), "_Clade"), "[", 2) # got this from: https://stackoverflow.com/questions/17215789/extract-a-substring-in-r-according-to-a-pattern
pullOut_noLF$ORspp <- sapply(strsplit(as.character(pullOut_noLF$Query), "_Clade"), "[", 1)
head(pullOut_noLF)
############# FIGURE OUT HITS THAT HAVE CHANGED ########
# NOTE THAT THIS INCLUDES THE 1-2 GENES THAT GOT SWITCHED DUE TO LEN 50 FILTER
# ARG.
missing <- pullOut_LengthFilter[!(pullOut_LengthFilter$bedname %in% pullOut_noLF$bedname),]
dim(missing) # 70. Hmm that's a lot. Only mising 1 functional.
pullOut_noLF$class <- "class" # standing in
dim(pullOut_LengthFilter)
dim(pullOut_noLF)
# for now: 
pullOut_all <- pullOut_noLF # this doesn't include the missing stuff! be ware!
############### GOING TO HAVE TO DEAL WITH THIS MORE FORMALLY I GUESS? #######

##################### FUNCTIONAL ORS ##############
functionalORs <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/mfur/step_9_FinalFunctionalList/mfur.classI.and.classII.final.ORs.passedInspection.txt",header=F) # this is the result of step 9 of functional analysis
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
exclude <- "mfur_GL897276.1_595480_596211_-1"
# and one functional gene isn't in pullOut
#missing
# find if any are missing in this not length filtered pullout:
missingFunc <- functionalORs_names[!(functionalORs_names %in% pullOut_all$bedname)]
# get the missing functional gene from the Functional length filter Pullout:

pullOut_allPlusMissingFunc <- rbind(pullOut_all[!(pullOut_all$bedname %in% exclude),],pullOut_LengthFilter[pullOut_LengthFilter$bedname==missingFunc,])
#####################################################

################# NON FUNCTIONAL ###############
################## a. TRUNCATED ###############
truncORs_NC <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/mfur/step_4_b_Visualize/VisualizationNotes/NC.missing.mfur.IDList.format.txt",header=F,strip.white=T) # this is list of truncated ORs. 
dim(truncORs_NC) #64
truncORs_N <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/mfur/step_4_b_Visualize/VisualizationNotes/Nmissing.Ccomplete.mfur.IDList.format.txt",header=F,strip.white=T)
dim(truncORs_N) #32
truncORs_C <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/mfur/step_4_b_Visualize/VisualizationNotes/Cmissing.Ncomplete.IDList.format.txt",header=F,strip.white=T)
dim(truncORs_C) # 11
# THIS IS NOT TRUNCATED: 
notTrunc_pgene <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/OR_PseudogeneResults/mfur/step_4_b_Visualize/VisualizationNotes/notTrunc.Pgene.IDList.format.txt",header=F,strip.white = T)

allTrunc <- unique(rbind(truncORs_C,truncORs_N,truncORs_NC)) # may not quite equal total because one may have doubled up if the ORFs were equally long -- using unique to eliminate 
# it's okay if sizes are close but don't quite match
sum(duplicated(rbind(truncORs_C,truncORs_N,truncORs_NC))) #1 that's where the dup is -- it's okay.  
colnames(allTrunc) <- "bedname"
dim(allTrunc)
################# b. possible pseudogenes (gotta dig into more) ################
# these are genes that are in pullOut_all, but aren't in truncated or in functional category. 
pgenes <- pullOut_allPlusMissingFunc[!(pullOut_allPlusMissingFunc$bedname %in% allTrunc$bedname) & !(pullOut_allPlusMissingFunc$bedname %in% functionalORs_names),]
dim(pgenes) # 377 
dim(allTrunc) # 106
length(unique(functionalORs_names)) # 816
# 377+106+816 = 1299
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
# maybe want to add a category -- there a bunch of pseudogenes that are really short
# but aren't called truncated bc not within 300bp of gap, or 30bp of ends of scaffold.
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$bedname %in% pgenes$bedname & pullOut_allPlusMissingFunc$Hit_length < 200,]$label <- "fragment"

# check:
sum(pullOut_allPlusMissingFunc$label=="pseudogene") # 377 (129 if fragments excluded) # why is this now 447, where is this extra 70 coming from?
sum(pullOut_allPlusMissingFunc$label=="functional") # (this should be 816 -- need to get MISSING back in here )
sum(pullOut_allPlusMissingFunc$label=="truncated") # 106 
# cool. 
  sum(pullOut_allPlusMissingFunc$label=="fragment") # 248 

#write.table(pullOut_allPlusMissingFunc,"step_5_labelAllgenes/pullOut.func.trunc.pgene.labelled.mfur.all.AddedMissing.fragments.txt",row.names=F,quote = F,sep="\t")
################# PLOT GENES BY CLADE ################
# # calculate % pseudogenized:
library(tidyr)
library(dplyr)
# fix awkward _CLASS_I
pullOut_allPlusMissingFunc$clade <- as.factor(pullOut_allPlusMissingFunc$clade)

levels(pullOut_allPlusMissingFunc$clade)[levels(pullOut_allPlusMissingFunc$clade)=="_CLASS_I"]<- "ClassI"
#pullOut_all[pullOut_all$clade=="_CLASS_I",]$clade <- "ClassI"
# https://stackoverflow.com/questions/40348854/geom-bar-percentage-for-each-group-in-x

tally1 <- pullOut_allPlusMissingFunc %>% group_by(clade,label) %>% tally %>% mutate(n/sum(n)) %>% mutate(sum(n))
# this is so cool! I need to spend more time with dplyr/tidyr
tally1
colnames(tally1) <- c("clade","label","n","perc","sumN")
head(tally1)
write.table(tally1,"step_5_labelAllgenes/tallyOfGenesByClade.category.txt",row.names = F,quote=F,sep="\t")

p2 <- ggplot(pullOut_allPlusMissingFunc,aes(x=clade,fill=label))+
  geom_bar(stat="count")+
  theme_bw()+
  scale_fill_manual(values=c(func,pg,tr))+
  ggtitle("Ferret Olfactory Receptors")+
  theme(legend.title = element_blank(),legend.position="top",legend.background = element_rect("transparent"))+
  geom_text(data=tally1[tally1$label=="pseudogene",], aes(x=clade,y=sumN+10,label=n),size=4,nudge_x = 0,nudge_y=0.6,color=pg)+
  theme(axis.text.x = element_text(angle=45, hjust=1)) 
p2

ggsave("step_5_labelAllgenes/labeledGeneBreakdown.countPgeneLabel.pdf",p2,device="pdf",width=7,height=4)

############ Class I and II ########################
pullOut_allPlusMissingFunc$class <- NA
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$clade=="ClassI",]$class <- "I"
pullOut_allPlusMissingFunc[pullOut_allPlusMissingFunc$clade!="ClassI",]$class <- "II"

## write out:
write.table(pullOut_allPlusMissingFunc,"step_5_labelAllgenes/pullOut.func.trunc.pgene.labelled.mfur.all.AddedMissing.Fragments.txt",row.names=F,quote = F,sep="\t")


tally2 <- pullOut_all %>% group_by(class,label) %>% tally %>% mutate(n/sum(n)) %>% mutate(sum(n))
# this is so cool! I need to spend more time with dplyr/tidyr
tally2
colnames(tally2) <- c("class","label","n","perc","sumN")
head(tally2)


