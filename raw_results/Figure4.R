library(reshape2)
library(hydroGOF)

#READ IN BULK DATA
###############################################################################
# FUNCTIONS USED BY get_statistics_for_ggplot

#Function which loads data, subsets it and orders it
data_processing<-function(path){
  results<-read.table(path)
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  results<-results[,colnames(results)=='ERR522956', drop=FALSE]
  return(results)
}

#Function to find correlation
correlation<-function(x,y) {
  (diag(cor(y,x,method="spearman")))
}

#Function that returns precision
make_precision<-function(ground_truth, tool_estimates, threshold_unexpr){
  TP<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates>threshold_unexpr])
  #print(TP)
  FP<-length(ground_truth[ground_truth<=threshold_unexpr & tool_estimates>threshold_unexpr])
  return(TP/(TP+FP))
}

#Function that returns precision value per cell
return_precision_per_cell<-function(truth_input_data, estimate_input_data){

  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_precision(truth_input_data[,j], estimate_input_data[,j], 0)
  }

  df<-do.call(cbind,results)
  colnames(df)<-colnames(truth_input_data)
  return(df)

}


#Function that returns recall
make_recall<-function(ground_truth, tool_estimates, threshold_unexpr){
  TP<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates>threshold_unexpr])
  #print(TP)
  FN<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates<=threshold_unexpr])
  return(TP/(TP+FN))
}

#Function that returns recall per cell
return_recall_per_cell<-function(truth_input_data, estimate_input_data){

  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_recall(truth_input_data[,j], estimate_input_data[,j], 0)
  }

  df<-do.call(cbind,results)
  colnames(df)<-colnames(truth_input_data)
  return(df)

}

#Function that returns F1
find_F1<-function(precision,recall){
  F1<-2*((precision*recall)/(precision + recall))
  return(F1)
}

#Function that returns legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}



################################################################################
#LOAD DATA AND PROCESS IT

#big function which loads data, processes it, finds statistics and returns it. Note this function assumes your results matrices are in a single directory and there are no other text files in that directory
get_statistics_for_ggplot<-function(path_to_dir,add_isoforms){
  #TODO: For given path, open + process all files in directory
  filenames <- list.files(path_to_dir, pattern="*_TPM.txt", full.names=TRUE)

  #create objects for files
  for (i in 1:length(filenames)){
    program<-strsplit(filenames[i], "/")[[1]][2]
    program<-strsplit(program,"clean_")[[1]][2]
    program<-strsplit(program,"_TPM.txt")[[1]][1]
    print(program)
    assign(program, data_processing(filenames[i]))
  }

  #If not RSEM, need to add unexpressed isoforms which were not input into the simulation process
  if (add_isoforms==TRUE){
    ground_truth<-read.table("data/BLUEPRINT_results_benchmark_2/polyester/results_matrices/clean_ground_truth_TPM.txt")
    ground_truth<-ground_truth[,colnames(ground_truth) %in% colnames(Salmon_align)]
    #return(ground_truth)
    extra_rows<-Salmon_align[!rownames(Salmon_align) %in% rownames(ground_truth),]
    extra_rows[extra_rows!=0]<-0
    ground_truth<-rbind(ground_truth,extra_rows)
    ground_truth<-ground_truth[ , order(colnames(ground_truth))]
    ground_truth<-ground_truth[order(rownames(ground_truth)),]
    #return(ground_truth)
  }

  #Find Spearman's rho for each method
  RSEM_cor<-correlation(RSEM,ground_truth)
  print(dim(ground_truth))
  print(dim(Salmon_align))
  Salmon_align_cor<-correlation(Salmon_align, ground_truth)
  Salmon_quasi_cor<-correlation(Salmon_quasi, ground_truth)
  Salmon_SMEM_cor<-correlation(Salmon_SMEM, ground_truth)
  Sailfish_cor<-correlation(Sailfish, ground_truth)
  eXpress_cor<-correlation(eXpress, ground_truth)
  Kallisto_cor<-correlation(Kallisto, ground_truth)

  #store in a dataframe
  spearmans_data<-melt(rbind(RSEM_cor,Salmon_align_cor, Salmon_quasi_cor, Salmon_SMEM_cor, Sailfish_cor, eXpress_cor, Kallisto_cor))
  #spearmans_data<-melt(rbind(Salmon_align_cor, Salmon_quasi_cor, Salmon_SMEM_cor, Sailfish_cor, eXpress_cor, Kallisto_cor))
  spearmans_data<-cbind("spearmans", spearmans_data)

  colnames(spearmans_data)<-c("Statistic","Tool", "Sample_name", "Value")

  RSEM_nmrse<-nrmse(log2(RSEM+1), log2(ground_truth +1))
  Salmon_align_nmrse<-nrmse(log2(Salmon_align+1), log2(ground_truth +1))
  Salmon_quasi_nmrse<-nrmse(log2(Salmon_quasi+1), log2(ground_truth +1))
  Salmon_SMEM_nmrse<-nrmse(log2(Salmon_SMEM+1), log2(ground_truth +1))
  Sailfish_nmrse<-nrmse(log2(Sailfish+1), log2(ground_truth +1))
  eXpress_nmrse<-nrmse(log2(eXpress+1), log2(ground_truth +1))
  Kallisto_nmrse<-nrmse(log2(Kallisto+1), log2(ground_truth +1))

  #Store in a dataframe
  nrmse_data<-melt(rbind(RSEM_nmrse,Salmon_align_nmrse, Salmon_quasi_nmrse, Salmon_SMEM_nmrse, Sailfish_nmrse, eXpress_nmrse, Kallisto_nmrse))
  nrmse_data<-cbind("nrmse", nrmse_data)
  colnames(nrmse_data)<-c("Statistic","Tool","Sample_name","Value")
  print(head(nrmse_data))

  #Find precision for each method
  RSEM_precision<-return_precision_per_cell(ground_truth,RSEM)
  Salmon_align_precision<-return_precision_per_cell(ground_truth,Salmon_align)
  Salmon_quasi_precision<-return_precision_per_cell(ground_truth,Salmon_quasi)
  Salmon_SMEM_precision<-return_precision_per_cell(ground_truth,Salmon_SMEM)
  Sailfish_precision<-return_precision_per_cell(ground_truth,Sailfish)
  eXpress_precision<-return_precision_per_cell(ground_truth,eXpress)
  Kallisto_precision<-return_precision_per_cell(ground_truth,Kallisto)

  #Store in dataframe
  precision_data<-rbind(RSEM_precision,Salmon_align_precision, Salmon_quasi_precision, Salmon_SMEM_precision, Sailfish_precision, eXpress_precision, Kallisto_precision)
  rownames(precision_data)<-c("RSEM_precision","Salmon_align_precision", "Salmon_quasi_precision", "Salmon_SMEM_precision", "Sailfish_precision", "eXpress_precision", "Kallisto_precision")
  precision_data<-melt(precision_data)
  precision_data<-cbind("precision", precision_data)
  colnames(precision_data)<-c("Statistic","Tool","Sample_name","Value")

  #Find recall for each method
  RSEM_recall<-return_recall_per_cell(ground_truth,RSEM)
  Salmon_align_recall<-return_recall_per_cell(ground_truth,Salmon_align)
  Salmon_quasi_recall<-return_recall_per_cell(ground_truth,Salmon_quasi)
  Salmon_SMEM_recall<-return_recall_per_cell(ground_truth,Salmon_SMEM)
  Sailfish_recall<-return_recall_per_cell(ground_truth,Sailfish)
  eXpress_recall<-return_recall_per_cell(ground_truth,eXpress)
  Kallisto_recall<-return_recall_per_cell(ground_truth,Kallisto)

  recall_data<-(rbind(RSEM_recall,Salmon_align_recall, Salmon_quasi_recall, Salmon_SMEM_recall, Sailfish_recall, eXpress_recall, Kallisto_recall))
  rownames(recall_data)<-c("RSEM_recall","Salmon_align_recall", "Salmon_quasi_recall", "Salmon_SMEM_recall", "Sailfish_recall", "eXpress_recall", "Kallisto_recall")
  recall_data<-melt(recall_data)
  recall_data<-cbind(statistic="recall", recall_data)
  colnames(recall_data)<-c("Statistic","Tool","Sample_name","Value")

  #Find F1 for each method
  RSEM_F1<-find_F1(RSEM_precision,RSEM_recall)
  Salmon_align_F1<-find_F1(Salmon_align_precision, Salmon_align_recall)
  Salmon_quasi_F1<-find_F1(Salmon_quasi_precision, Salmon_quasi_recall)
  Salmon_SMEM_F1<-find_F1(Salmon_SMEM_precision, Salmon_SMEM_recall)
  Sailfish_F1<-find_F1(Sailfish_precision, Sailfish_recall)
  eXpress_F1<-find_F1(eXpress_precision, eXpress_recall)
  Kallisto_F1<-find_F1(Kallisto_precision, Kallisto_recall)

  F1_data<-rbind(RSEM_F1, Salmon_align_F1, Salmon_quasi_F1, Salmon_SMEM_F1, Sailfish_F1, eXpress_F1, Kallisto_F1)
  rownames(F1_data)<-c("RSEM_F1", "Salmon_align_F1", "Salmon_quasi_F1", "Salmon_SMEM_F1", "Sailfish_F1", "eXpress_F1", "Kallisto_F1")
  F1_data<-melt(F1_data)
  F1_data<-cbind(statistic="F1", F1_data)
  colnames(F1_data)<-c("Statistic","Tool","Sample_name","Value")

  #Combine dfs into one
  df<-rbind(spearmans_data,nrmse_data,precision_data,recall_data, F1_data)
  return(df)
}

#DO STATISTICS
bulk_data<-get_statistics_for_ggplot("data", FALSE)

#SAVE DATA
write.table(bulk_data, "../figures/data/Figure4.txt")
