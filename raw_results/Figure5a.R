###############################################################################
# CREATE FILTERS BASED ON QC STATS

#QC script
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(gtable)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(scater)



#LOAD TEICHMAN DATA AND PROCESS IT
###############################################################################
# CREATE FILTERS BASED ON QC STATS

single_cell<-read.table("data/clean_Kallisto_real_Counts.txt")
single_cell<-single_cell[,colnames(single_cell)!='ERR522956', drop=FALSE]

rep_1<-read.table("data/filenames_2i_2.txt")
rep_2<-read.table("data/filenames_2i_3.txt")
rep_3<-read.table("data/filenames_2i_4.txt")
rep_4<-read.table("data/filenames_2i_5.txt")

#Function that formats batch data
make_batch<-function(counts_data){
  
  batch<-vector()
  
  batch_1<-colnames(counts_data) %in% rep_1$V1
  batch_2<-colnames(counts_data) %in% rep_2$V1
  batch_3<-colnames(counts_data) %in% rep_3$V1
  batch_4<-colnames(counts_data) %in% rep_4$V1
  
  for (i in 1:length(colnames(counts_data))){
    if (batch_1[i]==TRUE){
      batch[i]<-'2i_2'
    }else if (batch_2[i]==TRUE){
      batch[i]<-'2i_3'
    }
    else if (batch_3[i]==TRUE){
      batch[i]<-'2i_4'
    }
    else if (batch_4[i]==TRUE){
      batch[i]<-'2i_5'
    }
    else{
      print("Something went wrong")
    }
  }
  return(batch)
}

ids<-names(single_cell)
batch<-make_batch(single_cell)

anno<-as.data.frame(cbind(batch,ids))
rownames(anno)<-anno$ids

teich_scater <- SingleCellExperiment(
  assays = list(counts = as.matrix(single_cell)), 
  colData = anno
)

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
isSpike(teich_scater, "MT") <- rownames(teich_scater) %in% mt_isoforms

teich_scater_QC <- calculateQCMetrics(
  teich_scater,
  feature_controls = list(MT = isSpike(teich_scater, "MT"))
)

mt_reads<-plotPhenoData(
  teich_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_MT",
             colour = "batch")
)

teich_scater_QC<-teich_scater_QC[,teich_scater_QC$pct_counts_MT<10]

#Read in QC statistics files
QC_raw<-read.csv("data/raw/read_alignment_qc.csv", header=FALSE)
names(QC_raw)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")

#Filters applied based on Supplementary Figure 10
filter_ES<-QC_raw[QC_raw$NonUnique<2500000 & QC_raw$NumReads<12000000 & QC_raw$NumReads>3500000 & QC_raw$NumAlignments>4000000 & QC_raw$NumAlignments< 32000000,]
filter_ES<-filter_ES[filter_ES$Filename %in% teich_scater_QC$ids,]
filter_ES<-filter_ES[filter_ES$Filename!='ERR522956',]

rm(list=setdiff(ls(), c("filter_ES")))
#######################################

single_cell<-read.table("data/clean_ground_truth_counts.txt")

#filter_ES relevant cells
single_cell<-single_cell[,colnames(single_cell) %in% filter_ES$Filename]


ids<-names(single_cell)
batch<-rep("batch", ncol(single_cell))

anno<-as.data.frame(cbind(batch,ids))
rownames(anno)<-anno$ids

teich_scater <- SingleCellExperiment(
  assays = list(counts = as.matrix(single_cell)), 
  colData = anno
)

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
isSpike(teich_scater, "MT") <- rownames(teich_scater) %in% mt_isoforms

teich_scater_QC <- calculateQCMetrics(
  teich_scater,
  feature_controls = list(MT = isSpike(teich_scater, "MT"))
)

mt_reads<-plotPhenoData(
  teich_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_MT",
             colour = "batch")
)

teich_scater_QC<-teich_scater_QC[,teich_scater_QC$pct_counts_MT<10]

#Read in QC statistics files
QC_Sim<-read.csv("data/simulated/read_alignment_qc.csv")
names(QC_Sim)<-c("Filename","Unique","NonUnique","Unmapped","NumAlignments","NumReads")

#Keep only filtered cells
QC_Sim<-QC_Sim[QC_Sim$Filename %in% filter_ES$Filename,]

#Filters applied based on plots
filter_ES<-QC_Sim[QC_Sim$Unique>5000000,]
filter_ES<-filter_ES[filter_ES$Filename %in% teich_scater_QC$ids,]

rm(list=setdiff(ls(), c("filter_ES")))

data_processing<-function(path, filter_ES){
  results<-read.table(path)
  results<-results[,colnames(results) %in% filter_ES$Filename]
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  return(results)
}

#load data
ground_truth_ES<-data_processing("data/clean_ground_truth_TPM.txt", filter_ES)
RSEM_ES<-data_processing("data/clean_RSEM_TPM.txt", filter_ES)
Salmon_align_ES<-data_processing("data/clean_Salmon_align_TPM.txt", filter_ES)
Salmon_quasi_ES<-data_processing("data/clean_Salmon_quasi_TPM.txt", filter_ES)
Salmon_SMEM_ES<-data_processing("data/clean_Salmon_SMEM_TPM.txt", filter_ES)
Sailfish_ES<-data_processing("data/clean_Sailfish_TPM.txt", filter_ES)
eXpress_ES<-data_processing("data/clean_eXpress_TPM.txt", filter_ES)
Kallisto_ES<-data_processing("data/clean_Kallisto_TPM.txt", filter_ES)

#######################################################################################################

#Remove rows from ground truth matrix with more than 100%, 80%, 60%, 40%, 20% zeros
library(hydroGOF)

#Function to find correlation
correlation_mean<-function(x,y) {
  mean(diag(cor(y,x,method="spearman")))
}

#Function to find standard error in spearmans
correlation_std<-function(x,y) {
  std(diag(cor(y,x,method="spearman")))
}

#Function to find correlation
nrmse_mean<-function(estimates,truth) {
  mean(nrmse(log2(estimates+1),log2(truth+1),))
}

#Function to find standard error
std <- function(x) sd(x)/sqrt(length(x))

#Function to find standard error in spearmans
nrmse_std<-function(estimates,truth) {
  std(nrmse(log2(estimates+1),log2(truth+1)))
}

#Function to remove rows with more than number zeros
remove_zeros<-function(number,truth){
  return(truth[rowSums(truth==0)<(ncol(truth)*number),,drop=FALSE])
}

#Function to keep only rows in ground truth in expression estimates
filter<-function(estimate,truth) {
  return(estimate[rownames(estimate) %in% rownames(truth), , drop=FALSE])
}

#Function that returns misclassification rate
make_misclassify<-function(ground_truth, tool_estimates, threshold_unexpr){
  FN<-length(ground_truth[ground_truth>threshold_unexpr & tool_estimates<=threshold_unexpr])
  #print(TP)
  FP<-length(ground_truth[ground_truth<=threshold_unexpr & tool_estimates>threshold_unexpr])
  return((FN+FP)/length(ground_truth))
}

#Function that returns mean misclassification rate per cell
return_mean_misclassify_per_cell<-function(truth_input_data, estimate_input_data){

  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_misclassify(truth_input_data[,j], estimate_input_data[,j], 0)
  }

  return(mean(do.call(rbind,results)[,1]))

}

#Function that returns standard error of misclassification rate per cell
return_std_misclassify_per_cell<-function(truth_input_data, estimate_input_data){

  results<-list()
  for (i in 1:length(colnames(truth_input_data))){
    j=colnames(truth_input_data)[i]
    results[i]<-make_misclassify(truth_input_data[,j], estimate_input_data[,j], 0)
  }

  return(std(do.call(rbind,results)[,1]))

}

tools<-c("RSEM","Salmon Alignment", "Salmon Quasi", "Salmon SMEM", "Sailfish", "eXpress", "Kallisto")

#Function which takes a number, ground truth and expression estimates as input, removes rows with more than that percentage of zeros from ground truth, then returns stats
spearmans_rho<-function(number, cell_type, ground_truth,RSEM,Salmon_align, Salmon_quasi, Salmon_SMEM, Sailfish, eXpress, Kallisto){
  #Remove rows with more than number zeros
  ground_truth<-remove_zeros(number,ground_truth)

  #Keep only rows in ground truth in expression estimates
  RSEM<-filter(RSEM,ground_truth)
  Salmon_align<-filter(Salmon_align,ground_truth)
  Salmon_quasi<-filter(Salmon_quasi,ground_truth)
  Salmon_SMEM<-filter(Salmon_SMEM,ground_truth)
  Sailfish<-filter(Sailfish,ground_truth)
  eXpress<-filter(eXpress,ground_truth)
  Kallisto<-filter(Kallisto,ground_truth)

  #create empty vectors for results
  spearmans_results<- vector(mode="numeric", length=0)
  spearmans_error<- vector(mode="numeric", length=0)
  nrmse_results<-vector(mode="numeric", length=0)
  nrmse_error<-vector(mode="numeric", length=0)
  misclassify_results<-vector(mode="numeric", length=0)
  misclassify_error<-vector(mode="numeric", length=0)


  #Find spearmans rho
  spearmans_results[1]<-correlation_mean(RSEM,ground_truth)
  spearmans_results[2]<-correlation_mean(Salmon_align,ground_truth)
  spearmans_results[3]<-correlation_mean(Salmon_quasi,ground_truth)
  spearmans_results[4]<-correlation_mean(Salmon_SMEM, ground_truth)
  spearmans_results[5]<-correlation_mean(Sailfish,ground_truth)
  spearmans_results[6]<-correlation_mean(eXpress,ground_truth)
  spearmans_results[7]<-correlation_mean(Kallisto,ground_truth)

  #find standard error
  spearmans_error[1]<-correlation_std(RSEM,ground_truth)
  spearmans_error[2]<-correlation_std(Salmon_align,ground_truth)
  spearmans_error[3]<-correlation_std(Salmon_quasi,ground_truth)
  spearmans_error[4]<-correlation_std(Salmon_SMEM, ground_truth)
  spearmans_error[5]<-correlation_std(Sailfish,ground_truth)
  spearmans_error[6]<-correlation_std(eXpress,ground_truth)
  spearmans_error[7]<-correlation_std(Kallisto,ground_truth)

  #Find nrmse
  nrmse_results[1]<-nrmse_mean(RSEM,ground_truth)
  nrmse_results[2]<-nrmse_mean(Salmon_align,ground_truth)
  nrmse_results[3]<-nrmse_mean(Salmon_quasi,ground_truth)
  nrmse_results[4]<-nrmse_mean(Salmon_SMEM, ground_truth)
  nrmse_results[5]<-nrmse_mean(Sailfish,ground_truth)
  nrmse_results[6]<-nrmse_mean(eXpress,ground_truth)
  nrmse_results[7]<-nrmse_mean(Kallisto,ground_truth)

  #return vector of spearmans rho values
  nrmse_error[1]<-nrmse_std(RSEM,ground_truth)
  nrmse_error[2]<-nrmse_std(Salmon_align,ground_truth)
  nrmse_error[3]<-nrmse_std(Salmon_quasi,ground_truth)
  nrmse_error[4]<-nrmse_std(Salmon_SMEM, ground_truth)
  nrmse_error[5]<-nrmse_std(Sailfish,ground_truth)
  nrmse_error[6]<-nrmse_std(eXpress,ground_truth)
  nrmse_error[7]<-nrmse_std(Kallisto,ground_truth)

  return(cbind(tools,spearmans_results, spearmans_error, nrmse_results, nrmse_error, rep(number, 7), rep(cell_type, 7)))
}


hundred_ES<-spearmans_rho(1, "ES",ground_truth_ES, RSEM_ES,Salmon_align_ES, Salmon_quasi_ES, Salmon_SMEM_ES, Sailfish_ES, eXpress_ES, Kallisto_ES)
eighty_ES<-spearmans_rho(0.8,"ES",ground_truth_ES, RSEM_ES,Salmon_align_ES, Salmon_quasi_ES, Salmon_SMEM_ES, Sailfish_ES, eXpress_ES, Kallisto_ES)
sixty_ES<-spearmans_rho(0.6,"ES",ground_truth_ES, RSEM_ES,Salmon_align_ES, Salmon_quasi_ES, Salmon_SMEM_ES, Sailfish_ES, eXpress_ES, Kallisto_ES)
forty_ES<-spearmans_rho(0.4,"ES",ground_truth_ES, RSEM_ES,Salmon_align_ES, Salmon_quasi_ES, Salmon_SMEM_ES, Sailfish_ES, eXpress_ES, Kallisto_ES)
twenty_ES<-spearmans_rho(0.2,"ES",ground_truth_ES, RSEM_ES,Salmon_align_ES, Salmon_quasi_ES, Salmon_SMEM_ES, Sailfish_ES, eXpress_ES, Kallisto_ES)

figure_5a<-rbind(hundred_ES, eighty_ES, sixty_ES, forty_ES, twenty_ES)
colnames(figure_5a)[6]<-"percentage_zeros"
figure_5a[,6]<-as.numeric(figure_5a[,6]) * 100
colnames(figure_5a)[7]<-"cell_type"

write.table(figure_5a,"../figures/data/Figure5a.txt")
