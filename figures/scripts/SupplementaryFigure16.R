library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)

#######################################################################
# FUNCTIONS

#function removes trailing strings from tool names
remove_trail<-function(ggplot_results){
  ggplot_results$Tool<-sub("_cor", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_F1", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_nmrse", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_precision", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("_recall", "", ggplot_results$Tool)
  ggplot_results$Tool<-sub("Salmon_align", "Salmon\nAlign", ggplot_results$Tool)
  ggplot_results$Tool<-sub("Salmon_quasi", "Salmon\nQuasi", ggplot_results$Tool)
  ggplot_results$Tool<-sub("Salmon_SMEM", "Salmon\nSMEM", ggplot_results$Tool)
  return(ggplot_results)
}

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#function to create ggplot object
make_ggplot<-function(df, title, ylabel){
  p<-ggplot(df %>% dplyr::arrange(desc(Experiment)), aes(x=Tool, y=Value, colour=Experiment)) + geom_point( position=position_jitter(width = .2), stat = "identity") + facet_grid(~Tool, scales= "free_x",space = "free_x")
  p<-p + theme(legend.position = 'none', axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(), text = element_text(size=14), strip.text.x = element_text(size=14))
  p<-p + ggtitle(title) + ylab(ylabel) + scale_colour_manual(values=cbbPalette)
  return(p)
}

#################################

#TO DO

#READ IN BULK DATA
bulk<-read.table("../data/Figure4.txt")

#split into ES and Blueprint
ES_bulk<-bulk[bulk$Sample_name=="ERR522956",]

#delete Sample_name column and create an Experiment column
ES_bulk<-cbind(ES_bulk[,1:2], Value=ES_bulk[,4], Experiment="bulk")

#READ IN SINGLE CELL DATA
#Read in ES
ES_single<-read.table("../data/Figure2.txt")

#Sort out columns of ES_single
ES_single<-data.frame(Statistic=ES_single$statistic, Tool=ES_single$Var1, Value=as.numeric(ES_single$value), Experiment="single")

ES_df<-rbind(ES_bulk, ES_single)

#Remove trailing strings after tool names
ES_df<-remove_trail(ES_df)

#Figure 4 style ES plots
ES_spear<-ES_df[ES_df$Statistic=="spearmans",]
ES_nrmse<-ES_df[ES_df$Statistic=="nrmse",]
ES_precision<-ES_df[ES_df$Statistic=="precision",]
ES_recall<-ES_df[ES_df$Statistic=="recall",]
ES_F1<-ES_df[ES_df$Statistic=="F1",]


ES_spearmans<-make_ggplot(ES_spear, "Spearman's Rho", "Spearman's Rho")
ES_nrmse<-make_ggplot(ES_nrmse, "NRMSE", "NRMSE")
ES_precision<-make_ggplot(ES_precision, "Precision", "Precision")
ES_recall<-make_ggplot(ES_recall, "Recall", "Recall")
ES_F1<-make_ggplot(ES_F1, "F1", "F1")

ggarrange(ES_F1, ggarrange(ES_precision, ES_recall, nrow=2), ES_spearmans,ES_nrmse,                                         
          nrow = 2,
          ncol = 2,
          labels = c("A","","B", "C")                                 
) 

ggsave("../pdfs/SupplementaryFigure16.pdf", plot = last_plot(), width=170 *2, units=c("mm") )
ggsave("../pngs/SupplementaryFigure16.png", plot = last_plot(), width=170 *2, units=c("mm"))