library(ggplot2)
library(ggpubr)
library(reshape2)
library(hydroGOF)
library(scater)
library(tidyverse)
library(data.table)

###############################################################################
# QC FUNCTIONS


###############################################################################
# CREATE FILTERS BASED ON QC STATS

#QC script
library(ggplot2)
library(ggpubr)
library(reshape2)

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

#####################################################################
# BULK DATA PROCESSING

two_iso<-read.table("data/gene_isoform_relationship.txt")

get_two_iso<-function(path){
  #read in data
  df<-read.table(path,header=T)
  #keep only transcripts with 2 iso
  df<-df[df$target_id %in% two_iso$V2,]

  #Get a list of unique genes
  unique_genes<-unique(two_iso$V1)
  genes_to_keep<-list()
  j<-1

  for (i in 1:length(unique_genes)){
    #For each gene, find corresponding isoforms
    goi<-two_iso[two_iso$V1==unique_genes[i],]

    #find expression level of both isoforms
    expression_of_isoform_1<-df[df$target_id==as.character(goi$V2[1]),]$est_counts
    expression_of_isoform_2<-df[df$target_id==as.character(goi$V2[2]),]$est_counts

    #If both isoforms are expressed, add gene to genes_to_keep
    if (expression_of_isoform_1>2 & expression_of_isoform_2>2){
      genes_to_keep[j]<-as.character(unique_genes[i])
      j<-j+1
    }
  }

  #keep only genes in genes_to_keep in two_iso
  updated_two_iso<-two_iso[two_iso$V1 %in% genes_to_keep,]


  #keep only isoforms in updated_two_iso in df
  df<-df[df$target_id %in% updated_two_iso$V2,]

  #add gene name column to df
  #df<-df[order(df$target_id),]
  #updated_two_iso<-updated_two_iso[order(updated_two_iso$V2)]
  colnames(updated_two_iso)<-c("gene_name", "target_id")
  df<-merge(df,updated_two_iso)

  return(df)
}

bulk_ES<-get_two_iso("../Simulation/Kallisto_results_real_data/ERR522956/abundance.tsv")



###############################################################################
# Single cell processing

#Function which loads data, subsets it and orders it
data_processing<-function(path, filter,bulk){
  results<-read.table(path)
  results<-results[,colnames(results) %in% filter$Filename]
  colnames(results)<-sub("sample_", "Cell", colnames(results))
  colnames(results)<-sub("sample", "Cell", colnames(results))
  colnames(results)<-sub("bias_","", colnames(results))
  results<-results[ , order(colnames(results))]
  results<-results[order(rownames(results)),]
  return(results)
}

#Function which returns dataframe of how many isoforms are expressed
how_many_iso<-function(path_to_single_cell, filter,bulk,id) {

  #Read in single cell rna seq data
  single_cell<-data_processing(path_to_single_cell, filter)

  #Only keep isoforms in bulk
  single_cell<-single_cell[rownames(single_cell) %in% bulk$target_id,]


  genes<-unique(bulk$gene_name)
  no_expressed_iso<-data.frame(genes)


  #For each cell
  for (i in 1:ncol(single_cell)){

    cell_no_expressed_iso<-vector()

    #For each gene
    for (j in 1:length(genes)){
      isoform_1<-bulk[bulk$gene_name==as.character(genes[j]),]$target_id[1]
      isoform_1<-sum(single_cell[as.character(isoform_1), i]>2)

      isoform_2<-bulk[bulk$gene_name==as.character(genes[j]),]$target_id[2]
      isoform_2<-sum(single_cell[as.character(isoform_2), i]>2)

      cell_no_expressed_iso[j]<-isoform_1 + isoform_2

    }

    #Record number of expressed isoforms in dataframe
    no_expressed_iso<-cbind(no_expressed_iso, cell_no_expressed_iso)
    colnames(no_expressed_iso)[i+1]<-colnames(single_cell)[i]
    #return(names(single_cell[,i,drop=FALSE]))


  }


  #Return dataframe
  no_expressed_iso<-melt(no_expressed_iso)
  no_expressed_iso<-cbind(no_expressed_iso, ID=rep(id, nrow(no_expressed_iso)))
  return(no_expressed_iso)

}


ES_num_iso<-how_many_iso("data/clean_Kallisto_real_Counts.txt", filter_ES, bulk_ES, "ES")

write.table(ES_num_iso, "../figures/data/Figure6_number_of_isoforms.txt")


###########################################################################################

percentage<-function(numbers){
  return((sum(numbers)/length(numbers)) * 100)
}


percentage_vs_expression<-function(df,path,filter,bulk,id){
  
  #get percentage of genes expressing 2 isoforms
  no_zeros<-df #[df$value!=0,]
  no_zeros$value<-no_zeros$value==2
  two_iso_percentage<-no_zeros %>% group_by(genes) %>% summarise(percentage(value))

  #get expression of each gene
  single_cell<-data_processing(path,filter)
  #Only keep isoforms in bulk
  single_cell<-single_cell[rownames(single_cell) %in% bulk$target_id,]
  new_bulk<-data.frame(bulk$target_id, bulk$gene_name)
  new_bulk<-plyr::rename(new_bulk, c("bulk.target_id"= "rn", "bulk.gene_name"="genes"))
  setDT(single_cell, keep.rownames = TRUE)[]
  single_cell<-as.data.frame(single_cell)
  single_cell<-merge(single_cell, new_bulk)


  single_cell<-melt(single_cell)
  expression<-single_cell %>% group_by(genes) %>% summarise(sum(log2(value+1)))
  exprs_percent<-merge(expression,two_iso_percentage)
  exprs_percent<-cbind(exprs_percent, ID=rep(id, nrow(exprs_percent)))
  return(exprs_percent)
  #ggplot(data = exprs_percent, aes(y=`sum(log2(value + 1))`,x=`percentage(value)`)) + geom_point() + ylab("log2(Counts + 1)") + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")
}

ES_percent_exprs_plot<-percentage_vs_expression(ES_num_iso, "data/clean_Kallisto_real_Counts.txt", filter_ES, bulk_ES, "ES")


write.table(ES_percent_exprs_data, "../figures/data/Figure_6_percent_exprs.txt")
