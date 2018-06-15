###############################################################################
# CREATE FILTERS BASED ON QC STATS

#QC script
library(ggplot2)
library(ggpubr)
library(reshape2)
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





mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
ids<-names(single_cell)
batch<-make_batch(single_cell)

anno<-new("AnnotatedDataFrame", as.data.frame(cbind(batch,ids)))
rownames(anno)<-anno$ids
teich_scater <- scater::newSCESet(
  countData = single_cell,
  phenoData = anno
)

teich_scater_QC <- scater::calculateQCMetrics(
  teich_scater,
  feature_controls = list(MT = mt_isoforms)
)

mt_reads<-scater::plotPhenoData(
  teich_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_MT",
             colour = "batch")
)

teich_scater_QC<-teich_scater_QC[,teich_scater_QC$pct_counts_feature_controls_MT<10]

#Read in QC statistics files
QC_raw<-read.csv("../Simulation/QC_stats/raw/read_alignment_qc.csv")

#Filters applied based on plots
filter<-QC_raw[QC_raw$NonUnique<2500000 & QC_raw$NumReads<12000000 & QC_raw$NumReads>3500000 & QC_raw$NumAlignments>4000000 & QC_raw$NumAlignments< 32000000,]
filter<-filter[filter$Filename %in% teich_scater_QC$ids,]

rm(list=setdiff(ls(), c("filter")))
#######################################

single_cell<-read.table("data/clean_ground_truth_counts.txt")

#filter relevant cells
single_cell<-single_cell[,colnames(single_cell) %in% filter$Filename]

mt_isoforms<-c("ENSMUST00000082387", "ENSMUST00000082388", "ENSMUST00000082389", "ENSMUST00000082390", "ENSMUST00000082391", "ENSMUST00000082392", "ENSMUST00000082393", "ENSMUST00000082394", "ENSMUST00000082395", "ENSMUST00000082396", "ENSMUST00000082397", "ENSMUST00000082398", "ENSMUST00000082399", "ENSMUST00000082400", "ENSMUST00000082401", "ENSMUST00000082402", "ENSMUST00000082403", "ENSMUST00000082404", "ENSMUST00000082405", "ENSMUST00000082406", "ENSMUST00000082407", "ENSMUST00000082408", "ENSMUST00000082409", "ENSMUST00000082410", "ENSMUST00000082411", "ENSMUST00000082412", "ENSMUST00000084013", "ENSMUST00000082414", "ENSMUST00000082415", "ENSMUST00000082416", "ENSMUST00000082417", "ENSMUST00000082418", "ENSMUST00000082419", "ENSMUST00000082420", "ENSMUST00000082421", "ENSMUST00000082422", "ENSMUST00000082423")
ids<-names(single_cell)
batch<-rep("batch", ncol(single_cell))

anno<-new("AnnotatedDataFrame", as.data.frame(cbind(batch,ids)))
rownames(anno)<-anno$ids
teich_scater <- scater::newSCESet(
  countData = single_cell,
  phenoData = anno
)

teich_scater_QC <- scater::calculateQCMetrics(
  teich_scater,
  feature_controls = list(MT = mt_isoforms)
)

mt_reads<-scater::plotPhenoData(
  teich_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_feature_controls_MT",
             colour = "batch")
)

teich_scater_QC<-teich_scater_QC[,teich_scater_QC$pct_counts_feature_controls_MT<10]

#Read in QC statistics files
QC_Sim<-read.csv(""../Simulation/QC_stats/simulated/read_alignment_qc.csv")

#Keep only filtered cells
QC_Sim<-QC_Sim[QC_Sim$Filename %in% filter$Filename,]

#Filters applied based on plots
filter<-QC_Sim[QC_Sim$Unique>5000000,]
filter<-filter[filter$Filename %in% teich_scater_QC$ids,]

rm(list=setdiff(ls(), c("filter")))

################################################################################
#LOAD DATA AND PROCESS IT

#Function which loads data, subsets it and orders it
data_processing<-function(path, filter){
  results<-read.table(path)
  results<-results[,colnames(results) %in% filter$Filename]
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  return(results)
}



#load data
ground_truth<-data_processing("data/clean_ground_truth_TPM.txt", filter)

library(data.table)
setDT(ground_truth, keep.rownames = TRUE)[]


cells<-colnames(ground_truth)

ground_truth_expr<-ground_truth %>% gather(cells[2:ncol(ground_truth)], key="cell", value="estimates")
ground_truth_expr$estimates<-100*log2(ground_truth_expr$estimates +1)

write.table(ground_truth_expr, gzfile("../figures/data/SupplementaryFigure12.gz"))
