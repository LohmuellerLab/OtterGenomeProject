require(ggplot2)
require(scales)
require(RColorBrewer)
require(scales)
elutCol="#377EB8"
pbraCol="#4DAF4A"
elut2Col="#C77CFF"
################ read in results #############

results_elut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/MSMC_DemographicInference_WholeGenome/msmc/output_20180209_50_250DPFilter/sea_otter.msmc.out.final.txt",header=T)
results_pbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/MSMC_DemographicInference_WholeGenome/output_20180209_50_250DPFilter/giant_otter.msmc.out.final.txt",header=T)
results_elut2 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/MSMC_50_250DPFilter_simsInSppDIrs/output_20190215//elut2.msmc.out.final.txt",header = T) # this is just using lib 283, most up to date NSO results 20190216

################## set parameters ############
############## Mutation rate -- not adding NSO to this for now #################
grid_elut <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/genome-stats/mutationRate/elut.GridOfMutationRates.50_250DPFilter.10_14_20_MyaDivergence.txt",header=T)
grid_pbra <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/genome-stats/mutationRate/pbra.GridOfMutationRates.50_250DPFilter.10_14_20_MyaDivergence.txt",header=T)

grid_elut$spp <- "Sea Otter"
grid_pbra$spp <- "Giant Otter"

grids <- rbind(grid_elut,grid_pbra)

# choosing high, low, reasonable
# My choices: average between elut and pbra (they are very close)
highMu <- mean(grids[grids$divergence_yrs== 10400000 & grids$generation_time==10,]$mutationRate)
highMu
# High: 9mya + 10yr/gen (extreme high mutation rate) = average(3.558263e-08, 3.542660e-08)
# Reasonable:  11mya, 4yr/gen = average()
reasonableMu <- mean(grids[grids$divergence_yrs== 14400000 & grids$generation_time==4,]$mutationRate)
reasonableMu
# Low : 17 mya, 2yr/gen
lowMu <- mean(grids[grids$divergence_yrs== 20100000 & grids$generation_time==2,]$mutationRate)
lowMu
########################## generation time ###############
gen2=2
gen4=4
gen10=10


# Ralls uses 4, 8, 10 for elut gen times: http://www.sciencedirect.com/science/article/pii/000632078390037X/pdf?md5=6ad3b1bb4e7dac1d4a5b80c4649566be&pid=1-s2.0-000632078390037X-main.pdf
# Choices: 2 (goes with low mu), 4 (goes with reasonable mu), 8, 10 (goes with high mu)
######################## scale results #################
# to convert: from the msmc guide
# MSMC outputs times and rates scaled by the mutation rate per basepair per generation. First, scaled times are given in units of the per-generation mutation rate. This means that in order to convert scaled times to generations, divide them by the mutation rate. In humans, we used mu=1.25e-8 per basepair per generation.To convert generations into years, multiply by the generation time, for which we used 30 years.
# 
# To get population sizes out of coalescence rates, first take the inverse of the coalescence rate, scaledPopSize = 1 / lambda00. Then divide this scaled population size by 2*mu (yes, this factor 2 is different from the time scaling, sorry).

## make a function
# input is .final.txt file, read into a df
scaleMSMC <- function(input,mu,gen,spp,label,category="main"){
  input$Ne <- (1/input$lambda_00)/(2*mu) # note the factor of 2! (not in time scaling) confirmed correct: https://github.com/stschiff/msmc-tools/blob/master/plot_utils.py
  input$LeftYears <- gen*(input$left_time_boundary/mu)
  input$RightYears <- gen*(input$right_time_boundary/mu)
  input$label <- label
  input$spp <- spp
  input$category <- "main" # main result (not bootstrap)
  input$Left_generations <- (input$left_time_boundary/mu)
  input$Right_generations <- (input$right_time_boundary/mu)
  return(input)
}


################ Results scaled by 3 mutation rates (For Supplement) ############
results_elut_highMu <- scaleMSMC(results_elut,highMu,gen10,"S. Sea Otter",paste("S. Sea Otter\n3. high \u03BC = ",signif(highMu,3),", yr/gen = ",gen10,"\n",sep=""),"main")
results_elut_highMu$mu <- highMu

results_elut_reasonableMu <- scaleMSMC(results_elut,reasonableMu,gen4,"S. Sea Otter",paste("S. Sea Otter\n2. reasonable \u03BC = ",signif(reasonableMu,3),", yr/gen = ",gen4,"\n",sep=""),"main")
results_elut_reasonableMu$mu <- reasonableMu

results_elut_lowMu <- scaleMSMC(results_elut,lowMu,gen2,"S. Sea Otter",paste("S. Sea Otter\n1. low \u03BC = ",signif(lowMu,3),", yr/gen = ",gen2,"\n",sep=""),"main")
results_elut_lowMu$mu <- lowMu

results_pbra_highMu <- scaleMSMC(results_pbra,highMu,gen10,"Giant Otter",paste("Giant Otter\n3. high \u03BC  = ",signif(highMu,3),", yr/gen = ",gen10,"\n",sep=""),"main")
results_pbra_highMu$mu <- highMu

results_pbra_reasonableMu <- scaleMSMC(results_pbra,reasonableMu,gen4,"Giant Otter",paste("Giant Otter\n2. reasonable \u03BC = ",signif(reasonableMu,3),", yr/gen = ",gen4,"\n",sep=""),"main")
results_pbra_reasonableMu$mu <- reasonableMu

results_pbra_lowMu <- scaleMSMC(results_pbra,lowMu,gen2,"Giant Otter",paste("Giant Otter\n1. low \u03BC = ",signif(lowMu,3),", yr/gen = ",gen2,"\n",sep=""),"main")
results_pbra_lowMu$mu <- lowMu

# nso: 
results_elut2_highMu <- scaleMSMC(results_elut2,highMu,gen10,"N. Sea Otter",paste("N. Sea Otter\n3. high \u03BC = ",signif(highMu,3),", yr/gen = ",gen10,"\n",sep=""),"main")
results_elut2_highMu$mu <- highMu

results_elut2_reasonableMu <- scaleMSMC(results_elut2,reasonableMu,gen4,"N. Sea Otter",paste("N. Sea Otter\n2. reasonable \u03BC = ",signif(reasonableMu,3),", yr/gen = ",gen4,"\n",sep=""),"main")
results_elut2_reasonableMu$mu <- reasonableMu

results_elut2_lowMu <- scaleMSMC(results_elut2,lowMu,gen2,"N. Sea Otter",paste("N. Sea Otter\n1. low \u03BC = ",signif(lowMu,3),", yr/gen = ",gen2,"\n",sep=""),"main")
results_elut2_lowMu$mu <- lowMu

# put them all together for plotting:
results3Mu <- rbind(results_elut_highMu,results_elut_lowMu,results_elut_reasonableMu,results_pbra_highMu,results_pbra_lowMu,results_pbra_reasonableMu,results_elut2_highMu,results_elut2_lowMu,results_elut2_reasonableMu)

# put tables together for SI:
require(dplyr)
pbra_mer1 <- merge(results_pbra_highMu,results_pbra_reasonableMu, by=c("time_index","left_time_boundary","right_time_boundary","lambda_00","spp","category"), suffixes = c("",".reasonable"))
pbra_mer2 <- merge(pbra_mer1,results_pbra_lowMu,by=c("time_index","left_time_boundary","right_time_boundary","lambda_00","spp","category"),suffixes = c(".high",".low"))
#View(pbra_mer2)

write.table(arrange(subset(pbra_mer2,select = -c(category,label.low,label.high,label.reasonable,mu.low,mu.high,mu.reasonable,spp)),by=as.numeric(time_index)),"/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/MSMC_DemographicInference_WholeGenome/output_20180209_50_250DPFilter/pbra.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",quote=F,row.names=F,col.names=T,sep="\t")

elut_mer1 <- merge(results_elut_highMu,results_elut_reasonableMu, by=c("time_index","left_time_boundary","right_time_boundary","lambda_00","spp","category"), suffixes = c("",".reasonable"))
elut_mer2 <- merge(elut_mer1,results_elut_lowMu,by=c("time_index","left_time_boundary","right_time_boundary","lambda_00","spp","category"),suffixes = c(".high",".low"))
#View(elut_mer2)

write.table(arrange(subset(elut_mer2,select = -c(category,label.low,label.high,label.reasonable,mu.low,mu.high,mu.reasonable,spp)),by=as.numeric(time_index)),"/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/MSMC_DemographicInference_WholeGenome/output_20180209_50_250DPFilter/elut.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",quote=F,row.names=F,col.names=T,sep="\t")

# nso: 
elut2_mer1 <- merge(results_elut2_highMu,results_elut2_reasonableMu, by=c("time_index","left_time_boundary","right_time_boundary","lambda_00","spp","category"), suffixes = c("",".reasonable"))
elut2_mer2 <- merge(elut2_mer1,results_elut2_lowMu,by=c("time_index","left_time_boundary","right_time_boundary","lambda_00","spp","category"),suffixes = c(".high",".low"))
#View(elut_mer2)

write.table(arrange(subset(elut2_mer2,select = -c(category,label.low,label.high,label.reasonable,mu.low,mu.high,mu.reasonable,spp)),by=as.numeric(time_index)),"/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/MSMC_50_250DPFilter_simsInSppDIrs/output_20190215/elut2.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",quote=F,row.names=F,col.names=T,sep="\t")
################# SUPPLEMENTARY FIGURE : PLOT with 3 scalings (no bootstraps, makes plot too busy) #############
cairo_pdf("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/MSMC/3MutationRates/Elut.Pbra.plusNSO.High.Reasonable.Low.MSMC.20180209_50_250DPFilternewElut2Filter.10_14_20_MyaDivergence.Results.pdf",width = 8,height=5)
p3Mu <- ggplot(results3Mu,aes(x=LeftYears,y=Ne,color=label,group=label))+
  geom_step(stat="identity",size=1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(paste("Sea Otter & Giant Otter PSMC' Results\nDifferent Mutation Rate & Generation Time Scalings",sep=""))+
  xlab("Years Ago")+
  ylab("Inverse Coalescence Rate, scaled by 2\u03BC\n(Ne in panmictic population)")+
  scale_y_log10(breaks=c(100,1000,10000,100000,1000000),labels=comma)+
  scale_x_log10(breaks=c(1000,2000,5000,10000,seq(25000,50000,by=25000),100000,200000,500000,1000000),labels=comma)+
  theme(legend.position= c(0.38,0.8),legend.background = element_rect(fill="transparent"))+
  scale_color_manual(values=c(c(brewer.pal(4,"Greens")[2:4]),c(brewer.pal(4,"Purples")[2:4]),c(brewer.pal(4,"Blues")[2:4])))+
  theme(legend.direction=("vertical"),legend.position=c(0.33,0.77),legend.background = element_rect(fill="transparent"),legend.text=element_text(size=8),legend.key.size = unit(.28,"cm"))+
  guides(colour = guide_legend(ncol=2,byrow=FALSE))
p3Mu
dev.off()


############################# For main text: just reasonable mu ###############
############### Bootstraps (just scale by reasonable mu) ############
numBoot=20


############# PBRA bootstraps ###############

pbra.path="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGiantOtterReadsToFerret/MSMC_DemographicInference_WholeGenome/output_20180209_50_250DPFilter/bootstrap-output" # path to files
pbra.boots <- list.files(path=pbra.path,pattern = "final.txt") # list of .final.txt files
# read in the bootstraps:
for (i in 1:(length(pbra.boots))){
  # import the file
  pbra.file <- read.table(file = paste(pbra.path,"/",pbra.boots[i],sep=""),header=T)
  pbra.file <- scaleMSMC(pbra.file,reasonableMu,gen=4,"Giant Otter","bootstrap")
  pbra.file$bootNum <- i
  my.name <- paste("pbra.boot_",i,sep="")
  # assign the name to the object
  assign(paste(my.name), pbra.file) # this is awesome!! makes each one its own df
  # scale:
}
# there are now 20 dfs named:
# pbra.boot_1 , pbra.boot_2 etc. ... pbra.boot_20
#combine them:
bootstrap_allPbra_reasonable4 <- pbra.boot_1
for(i in seq(2,20)){
  my.name <- paste("pbra.boot_",i,sep="")
  bootstrap_allPbra_reasonable4 <- rbind(bootstrap_allPbra_reasonable4,get(my.name))
}


############# ELUT bootstraps ###############
# 20180212: updated to be 50/250 DP bootstraps

elut.path="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/MSMC_DemographicInference_WholeGenome/msmc/output_20180209_50_250DPFilter/bootstrap-output/" # path to files
elut.boots <- list.files(path=elut.path,pattern = "final.txt") # list of .final.txt files
for (i in 1:(length(elut.boots))){
  # import the file
  elut.file <- read.table(file = paste(elut.path,"/",elut.boots[i],sep=""),header=T)
  elut.file <- scaleMSMC(elut.file,reasonableMu,gen=4,"S. Sea Otter","bootstrap")
  elut.file$bootNum <- i
  my.name <- paste("elut.boot_",i,sep="")
  # assign the name to the object
  assign(paste(my.name), elut.file) # this is awesome!! makes each one its own df
  # scale:
}
# there are now 20 dfs named:
# elut.boot_1 , elut.boot_2 etc. ... elut.boot_20
#combine them:
bootstrap_allElut_reasonable4 <- elut.boot_1
for(i in seq(2,20)){
  my.name <- paste("elut.boot_",i,sep="")
  bootstrap_allElut_reasonable4 <- rbind(bootstrap_allElut_reasonable4,get(my.name))
}

# nso:
elut2.path="/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/MSMC_50_250DPFilter_simsInSppDIrs/output_20190215/bootstrap-output/" # path to files
elut2.boots <- list.files(path=elut2.path,pattern = "final.txt") # list of .final.txt files
for (i in 1:(length(elut2.boots))){
  # import the file
  elut2.file <- read.table(file = paste(elut2.path,"/",elut2.boots[i],sep=""),header=T)
  elut2.file <- scaleMSMC(elut2.file,reasonableMu,gen=4,"N. Sea Otter","bootstrap")
  elut2.file$bootNum <- i
  my.name <- paste("elut2.boot_",i,sep="")
  # assign the name to the object
  assign(paste(my.name), elut2.file) # this is awesome!! makes each one its own df
  # scale:
}
# there are now 20 dfs named:
# elut2.boot_1 , elut2.boot_2 etc. ... elut2.boot_20
#combine them:
bootstrap_allelut2_reasonable4 <- elut2.boot_1
for(i in seq(2,20)){
  my.name <- paste("elut2.boot_",i,sep="")
  bootstrap_allelut2_reasonable4 <- rbind(bootstrap_allelut2_reasonable4,get(my.name))
}
################# combine species #################
results_reasonable4 <- rbind(results_elut_reasonableMu,results_pbra_reasonableMu,results_elut2_reasonableMu)
bootstrapResults_reasonable4 <- rbind(bootstrap_allElut_reasonable4,bootstrap_allPbra_reasonable4,bootstrap_allelut2_reasonable4)
################ PLOTS ###################

################# MAIN TEXT FIGURE ::: PLOT Inv Coal Rate and generations (reasonable mu), WITHTRIMMING ##########
# Trimming information: Chose index 33 for elut, 20 for pbra:
elutTimeIndex=33
pbraTimeIndex=20
elut2TimeIndex=31 # not the same as elut1 after new filtering

# because you use left generatins and are setting that time index you remove everything before
# the time index -1 time index. (based on python script)
# so new Na is:
############# NOTE index -1: (spent a lot of time tirekicking; this is right)
results_reasonable4[results_reasonable4$time_index==(elutTimeIndex-1) & results_reasonable4$spp=="S. Sea Otter",c("Left_generations","Ne")] # this is if you use elutTimeIndex in the python script; you trim off things before timeIndex-1 (because using left_generations)
results_reasonable4[results_reasonable4$time_index==(pbraTimeIndex-1) & results_reasonable4$spp=="Giant Otter",c("Left_generations","Ne")]
results_reasonable4[results_reasonable4$time_index==(elut2TimeIndex-1) & results_reasonable4$spp=="N. Sea Otter",c("Left_generations","Ne")]

############ NOTE A CHANGE IN DEPICTING CUTOFFS: need to draw a line at (index -1) because you are totally removing that index in python

# re order levels so key goes go --> sso --> nso to be consistent with other figures
levels(as.factor(results_reasonable4$spp))
results_reasonable4$spp <- factor(results_reasonable4$spp,levels=c("Giant Otter" , "S. Sea Otter" ,"N. Sea Otter"))
pFinal <- ggplot(results_reasonable4,aes(x=Left_generations,y=Ne,color=spp))+
  geom_rect(aes(xmin=results_reasonable4[results_reasonable4$spp=="S. Sea Otter" & results_reasonable4$time_index==elutTimeIndex-1,]$Left_generations, xmax=Inf, ymin=0, ymax=Inf),fill="gray80",alpha=0.01,color="transparent")+
  geom_step(stat="identity",size=1)+
  geom_step(stat="identity",size=0.5,alpha=0.2,data=bootstrapResults_reasonable4,aes(x=Left_generations,y=Ne,group=interaction(spp,bootNum),color=spp))+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle(paste("S. Sea Otter, N. Sea Otter & Giant Otter PSMC' Results\n\u03BC = ",signif(reasonableMu,3)," mutations/bp/gen",sep=""))+
  xlab("Generations Ago")+
  ylab("Inverse Coalescence Rate, scaled by 2\u03BC\n(Ne in panmictic population) ")+
  scale_y_log10(labels=comma)+
  scale_x_log10(labels=comma)+
  theme(legend.position= c(0.6,0.85),legend.background = element_rect(fill="transparent"))+
  scale_color_manual(values=c(pbraCol,elutCol,elut2Col))+
  geom_vline(xintercept = results_reasonable4[results_reasonable4$spp=="S. Sea Otter" & results_reasonable4$time_index==elutTimeIndex-1,]$Left_generations,linetype=2,color="grey20")+
  theme(legend.direction=("vertical"),legend.position=c(0.5,0.8),legend.background = element_rect(fill="transparent"),legend.text=element_text(size=14),legend.key.size = unit(1,"cm"))
pFinal


ggsave("/Users/annabelbeichman/Documents/UCLA/Otters/GidgetDeNovoGenome/Genome_Manuscript/Figures/MSMC/MainMSMCPlot/Elut.Pbra.plusNSO.Ne_IICR.Generations.Trim.TimeIndex34.20.ProperLineDrawn.50_250DPFilter.newElut2.Filter.14_MyaDivergence.pdf",pFinal,width = 6.75,height=5,device=cairo_pdf)



  

  
