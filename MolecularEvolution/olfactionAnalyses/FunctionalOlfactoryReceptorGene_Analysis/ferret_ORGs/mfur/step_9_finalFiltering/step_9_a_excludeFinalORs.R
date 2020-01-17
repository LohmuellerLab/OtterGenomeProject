# after step 8's visual inspection, have a list to exclude
setwd("/Users/annabelbeichman/Documents/UCLA/Otters/Olfaction/20170707_MainOlfactionResults/mfur/")
classI <- read.table("step_8_alignClassesSeparately/classI.tips.sppOnly.noOutlier.txt",header=F)
classII <- read.table("step_8_alignClassesSeparately/classII.tips.sppOnly.noOutlier.txt",header=F)

exclude <- read.table("step_8_alignClassesSeparately/step_8_Exclude.txt",header=F)
### make sure you removed the / after the slash stuff at the end of the file!!! 
# otherwise this wont work!!
classI.final <- classI[!(classI$V1 %in% exclude$V1),]
classII.final <- classII[!(classII$V1 %in% exclude$V1),]
dim(classI)
length(classI.final) # make sure lengths changed!!
dim(classII)
length(classII.final)
dim(classII)


### write them out as separate classes and together:
 write.table(classI.final,"step_9_FinalFunctionalList/classI.final.ORs.passedInspection.txt",quote=F,row.names=F,col.names=F)
 write.table(classII.final,"step_9_FinalFunctionalList/classII.final.ORs.passedInspection.txt",quote=F,row.names=F,col.names=F)
write.table(c(as.character(classI.final),as.character(classII.final)),"step_9_FinalFunctionalList/classI.and.classII.final.ORs.passedInspection.txt",quote=F,row.names=F,col.names=F)
