# gather up significant results to visualize 
gatherSignificantResults_scriptMaker <- function(sigResults,branch,date,rep){
  inspect <- sigResults # whatever subset of the genes you want to inspect 
  # fix group names and cluster names by adding leading zeros
  inspect$group <- sprintf("%02d",inspect$group)
  inspect$cluster <- sprintf("%05s",inspect$cluster)
  # get rid of slashes 
  inspect$Symbol <- gsub("N/A","N_A",inspect$Symbol)
  # set up directories
  write(paste("mkdir -p $SCRATCH/visualInspection/",branch,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"_",inspect$Symbol,sep=""),file=paste("output_codeml_",date,"/inspect_codeml_subset_",date,"_",branch,"_rep",rep,".sh",sep=""),append=T)
  # cp beb info 
  write(paste("cp -r group_",inspect$group,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"/codeml_branchSite_results_rep_",rep,"_",date,"/bebInfo $SCRATCH/visualInspection/",branch,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"_",inspect$Symbol,sep=""),file=paste("output_codeml_",date,"/inspect_codeml_subset_",date,"_",branch,"_rep",rep,".sh",sep=""),append=T)
  # cp original guidance html and gblocks html
  write(paste("cp -r group_",inspect$group,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"/guidance*/MSA*htm* $SCRATCH/visualInspection/",branch,"/cluster_",inspect$cluster,"_group_",inspect$group,"_",inspect$TName,"_",inspect$Symbol,sep=""),file=paste("output_codeml_",date,"/inspect_codeml_subset_",date,"_",branch,"_rep",rep,".sh",sep=""),append=T)
}

inspect_3 <- otters_3[otters_3$TName %in% int231,]
gatherSignificantResults_scriptMaker(sigResults = inspect_3,branch="otters",date = 20170705,rep=1)

