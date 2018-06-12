library(ggplot2)
library(ggpubr)
library(reshape2)
library(hydroGOF)
library(scater)

####################################
# FUNCTIONS

#Function for making graphs of results
plot_data<-function( df, title, ylabel, xlabel) {
  df$Var1 <- as.character(df$Var1)
  df$Var1 <- factor(df$Var1, levels=unique(df$Var1))
  
  spearmans<-ggplot(data=df, aes(x=Var1, y=value)) + geom_jitter(alpha=0.5, position=position_jitter(width = .2), aes(colour=Var1))  + stat_summary(fun.y=mean, geom="point", shape=95, size = 20, colour="black") 
  spearmans<- spearmans + scale_x_discrete(labels=c("RSEM", "Salmon Alignment", "Salmon Quasi", "Salmon SMEM", "Sailfish", "eXpress", "Kallisto"))
  spearmans<-spearmans + theme(axis.text.x=element_text( angle=30,vjust=.8, hjust=0.8), legend.position = 'none', text = element_text(size=14)) + scale_colour_manual(values=cbbPalette) + ylab(ylabel) + ggtitle(title) + xlab(xlabel)
  return(spearmans)
}

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##############################################
# MAKE FIGURE 2

#Read in performance statistics
ggplot_results<-read.table("../data/Figure2.txt")

#Split ggplot_results by statistic
spearmans<-ggplot_results[ggplot_results$statistic=="spearmans",]
nrmse<-ggplot_results[ggplot_results$statistic=="nrmse",]
precision<-ggplot_results[ggplot_results$statistic=="precision",]
recall<-ggplot_results[ggplot_results$statistic=="recall",]
F1<-ggplot_results[ggplot_results$statistic=="F1",]

#make graphs
spearmans_graph<-plot_data(spearmans, "Spearman's Rho", "Spearman's Rho", "")
nrmse_graph<-plot_data(nrmse, "NRMSE", "NRMSE", "")
precision_graph<-plot_data(precision,"Precision", "Precision", "")
recall_graph<-plot_data(recall, "Recall", "Recall", "")
F1_graph<-plot_data(F1, "F1", "F1", "")

ggarrange(F1_graph, ggarrange(precision_graph, recall_graph, nrow=2), spearmans_graph,nrmse_graph,                                         
          nrow = 2,
          ncol = 2,
          labels = c("A","","B", "C")                                 
) 
ggsave("../pdfs/Figure2.pdf", plot=last_plot(), height= 225, width=170, units=c("mm"))
ggsave("../pngs/Figure2.png", plot=last_plot(), height= 225, width=170, units=c("mm"))