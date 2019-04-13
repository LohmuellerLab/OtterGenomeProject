
# DP dist for top 50 scaffs: this is the result of steps 10a-b-c (DP dist per first 50 scaffs, then bin it per scaffold just to make file smaller)
dp50 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/stats/50Scaffs_DP_SumsPerScaff.txt",header=F)
colnames(dp50) <- c("group","sumDP","sites")

meanDP <- sum(as.numeric(dp50$sumDP))/sum(as.numeric(dp50$sites))
meanDP # 56 for sea otter 1 (gidget); 34 for giant otter;  39 for northern sea otter (was 73 with both libraries)

minDP <- meanDP*0.5 # 50% of mean DP
maxDP <- meanDP * 2.5 # 250% of mean DP ** be sure to not use 150% ** 
minDP
maxDP
cat(paste("meanDP: ",meanDP,"\nminDP: ",minDP,"\nmaxDP: ",maxDP,sep=""),file="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/stats/minMaxDP.elut2.283only.txt")

# plot DP distribution:
scaff1 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/stats/DP.dist.scaff.1.out.gz",header=F)

p1 <- ggplot(scaff1,aes(x=V1))+
  geom_density()+
  scale_x_continuous(limits=c(0,200))+
  geom_vline(xintercept = minDP,linetype="dashed",color="red")+
  geom_vline(xintercept = maxDP,linetype="dashed",color="red")+
  theme_bw()+
  xlab("DP")+
  ggtitle("Scaffold 1 DP Distribution")
p1
ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/stats/DP.dist.scaff1.png",device="png",p1,height=5,width=7)
