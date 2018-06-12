#QC script
library(ggplot2)
library(ggpubr)
library(scater)

single_cell<-read.table("Dropbox/Teichmann/clean_Kallisto_real_Counts.txt")
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
save(teich_scater_QC, file="../figures/data/SupplementaryFigure10_scater_object.RData", compress="bzip2", compression_level = 9)

#Read in QC statistics files
QC_raw<-read.csv("data/Teichmann_results-master/raw/read_alignment_qc.csv")
QC_raw<-QC_raw[QC_raw$Filename!='ERR522956',]

write.table(QC_raw, "../figures/data/SupplementaryFigure10_reads_alignment_data.txt")
