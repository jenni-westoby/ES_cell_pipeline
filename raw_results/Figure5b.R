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

###########################################################################################################


plot_cor_expr_zeros<-function(path,filter_filter, x_lab_text, y_lab_text, title, legend, cell_names){
  ground_truth<-data_processing(path,filter_filter)
  percent_zeros<-apply(ground_truth, 1, FUN=function(x) (length(x[x==0])/length(x))*100)
  expression<-rowMeans(log2(ground_truth+1))
  correlation_zeros_exprs<-data.frame(percent_zeros, expression)
  correlation_zeros_exprs<-data.frame(correlation_zeros_exprs, rep(cell_names, nrow(correlation_zeros_exprs)))
  return(correlation_zeros_exprs)
}

ES_sim_cor_expr_zeros<-plot_cor_expr_zeros("data/clean_ground_truth_counts.txt",filter_ES," ", " ", "\nSimulated ES cells", FALSE, "simulated_ES")
ES_real_cor_expr_zeros<-plot_cor_expr_zeros("data/clean_Kallisto_real_Counts.txt", filter_ES, "log2(counts + 1)", " ", "\nReal ES cells", FALSE, "real_ES")

df<-rbind(ES_sim_cor_expr_zeros, ES_real_cor_expr_zeros)
colnames(df)[3]<-"ID"
write.table(df, "../figures/data/Figure5b.txt")
#######################################################################################################
