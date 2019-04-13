### Making heterozygosity table (ignore Excel ones because they are potentially messedup)
require(ggplot2)
require(ggpubr)
jarOriginal <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/HeterozygosityTables/RobinsonHetTable.DrawFrom.txt",header=T,quote="",sep="\t",strip.white = T)
# had to save as a mac os formatted text file
head(jarOriginal)

# Restrict to genome-wide and animal only:
jarOriginal_gw_animal <- jarOriginal[jarOriginal$Seq_Type=="genome-wide" & jarOriginal$Kingdom=="animal",]
dim(jarOriginal)
dim(jarOriginal_gw_animal)
# just have one entry per common name:
jarOriginal_gw_animal_unique <- jarOriginal_gw_animal[!duplicated(jarOriginal_gw_animal$Common_Name),] # 68 ; one entry per species
dim(jarOriginal_gw_animal_unique)
# merge in IUCN rankings:
IUCN <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/HeterozygosityTables/IUCN_Categories.txt",header=T,sep="\t",quote="",strip.white = T)
dim(IUCN)

# merge IUCN with jar:
hetWIUCN <- merge(IUCN, jarOriginal_gw_animal_unique,by="Common_Name",all.x=F,all.y=F)
dim(hetWIUCN)
# get rid of duplicate entries **BE VERY CAREFUL. BAR PLOT WITH MULTIPLE ENTRIES WILL ADD TOO MUCH TO YOUR BARS (WILL ADD ALL TOGETHER FOR THOSE ENTRIES.) *MAKE SURE THIS ISN'T AN ISSUE ELSEWHERE

hetWIUCN <- unique(hetWIUCN)
dim(hetWIUCN)


# Make some additions: island fox (san nic), gray fox, sea otter, giant otter

additions <- read.table("~/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/HeterozygosityTables/AdditionsToHetTable.50_250DPFilter.txt",header=T,sep="\t")
# add to het:
colnames(additions) <- c("Common_Name","Species.y","IUCN_Status","Kingdom","Observed_pi","Seq_Type","Source","DOI","Plot_Name","Category")
hetWIUCN_plus <- rbind(additions,hetWIUCN[,c("Common_Name","Species.y","IUCN_Status","Plot_Name","Category","Kingdom","Observed_pi","Seq_Type","Source","DOI")])
# make some prelim plots:

# reorder labels:
require(dplyr)
levels(hetWIUCN_plus$IUCN_Status)
# [1] "CR" "EN" "LC" "NE" "NT" "VU"
hetWIUCN_reorder <- hetWIUCN_plus %>% mutate(IUCN_Status = factor(IUCN_Status, levels = c("NE","LC","NT","VU","EN","CR")))

require(ggplot2)
require(RColorBrewer)


#### WRITE OUT version of table for suppplement #### 
hetWIUCN_forSupplement <- hetWIUCN_reorder[hetWIUCN_reorder$Category=="mammal",c("Common_Name", "Species.y" ,  "IUCN_Status" , "Observed_pi" , "Seq_Type", "Source"    ,  "DOI" )]
colnames(hetWIUCN_forSupplement) <- c("Common_Name", "Species" ,  "IUCN_Status" , "Observed_pi" , "Seq_Type", "Source"    ,  "DOI" )
write.table(hetWIUCN_forSupplement,"/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Tables/HeterozygosityTables/HetTableForSupplement_matchesFigure.50_250DPFilter.justMammals.txt",sep="\t",row.names=F,quote=F)

########## Make a lollipop plot #########
hetWIUCN_reorder$IUCN_Long <- "NA"
hetWIUCN_reorder[hetWIUCN_reorder$IUCN_Status=="VU",]$IUCN_Long <- "Vulnerable (VU)"
hetWIUCN_reorder[hetWIUCN_reorder$IUCN_Status=="LC",]$IUCN_Long <- "Least Concern (LC)"
hetWIUCN_reorder[hetWIUCN_reorder$IUCN_Status=="EN",]$IUCN_Long <- "Endangered (EN)"
hetWIUCN_reorder[hetWIUCN_reorder$IUCN_Status=="NE",]$IUCN_Long <- "Not Evaluated (NE)"
hetWIUCN_reorder[hetWIUCN_reorder$IUCN_Status=="CR",]$IUCN_Long <- "Critcally Endangered (CR)"
hetWIUCN_reorder[hetWIUCN_reorder$IUCN_Status=="NT",]$IUCN_Long <- "Near Threatened (NT)"

require(dplyr)
hetWIUCN_reorder2 <- hetWIUCN_reorder %>% mutate(IUCN_Long = factor(IUCN_Long, levels = c("Not Evaluated (NE)","Least Concern (LC)","Near Threatened (NT)","Vulnerable (VU)","Endangered (EN)","Critcally Endangered (CR)")))

p3 <- ggdotchart(hetWIUCN_reorder2[hetWIUCN_reorder2$Category=="mammal",],x="Plot_Name",y="Observed_pi",color="IUCN_Long", palette=c("gray",brewer.pal(name="Oranges",5)),ggtheme = theme_bw(),add="segments",dot.size=10,label="IUCN_Status",font.label = list(color = "black", size = 14, vjust=0.5), add.params = list(color = "lightgray", size = 2))+
  coord_flip()+
  ylab(expression(paste("Observed ",pi)))+
  xlab("Species")+
  theme(axis.text.x=element_text(size=18),axis.title=element_text(size=18),axis.text.y=element_text(size=18),title=element_text(size=18))+
  theme(legend.position = c(0.8, 0.5),legend.background = element_rect("transparent"),legend.text=element_text(size=16),legend.title=element_blank())+
  ggtitle("Genome-Wide Heterozygosity in Mammals")
p3
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/Fig1.GenomeWideHeterozygosity.justMammals.50_250DPFilter.LOLLIPOP.inclNSO.pdf",p3,device = "pdf",height = 10,width = 14,units = "in")

##################### request from MBE : exclude house mouse ###########
hetWIUCN_reorder3 <- hetWIUCN_reorder2[hetWIUCN_reorder2$Plot_Name!="House mouse",]
p4 <- ggdotchart(hetWIUCN_reorder3[hetWIUCN_reorder3$Category=="mammal",],x="Plot_Name",y="Observed_pi",color="IUCN_Long", palette=c("gray",brewer.pal(name="Oranges",5)),ggtheme = theme_bw(),add="segments",dot.size=10,label="IUCN_Status",font.label = list(color = "black", size = 14, vjust=0.5), add.params = list(color = "lightgray", size = 2))+
  coord_flip()+
  ylab(expression(paste("Observed ",pi)))+
  xlab("Species")+
  theme(axis.text.x=element_text(size=18,angle=0,hjust=0.5),axis.title=element_text(size=18),axis.text.y=element_text(size=18),title=element_text(size=18))+
  theme(legend.position = c(0.8, 0.5),legend.background = element_rect("transparent"),legend.text=element_text(size=16),legend.title=element_blank())+
  ggtitle("Genome-Wide Heterozygosity in Mammals")
p4
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/Heterozygosity/Fig1.GenomeWideHeterozygosity.justMammals.50_250DPFilter.LOLLIPOP.NoHouseMouse.inclNSO.USETHIS.pdf",p4,device = "pdf",height = 11,width = 11,units = "in")
