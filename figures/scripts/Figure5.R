library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(gtable)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(MASS)
library(viridis)

####################################
# FUNCTIONS

#Function to process data
process_data<-function(df,cell_type){
  df<-df[df$cell_type==as.name(cell_type),]
  df<-cbind(df, max_spearmans_error = df$spearmans_results + df$spearmans_error, min_spearmans_error = df$spearmans_results - df$spearmans_error)
  df<-cbind(df, max_nrmse_error = df$nrmse_results + df$nrmse_error, min_nrmse_error = df$nrmse_results - df$nrmse_error)
  return(df)
}

#Function to make figure legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

#Function to plot results
plot_results<-function(title, yaxis_text, xaxis_text, legend_true, df, results,max_error, min_error){

  #Plot results
  p1<-ggplot(data=df, aes_(x= ~percentage_zeros, y=as.name(results), group=~tools, colour=~factor(tools)))+ geom_line() + geom_errorbar(aes_(ymin=as.name(min_error), ymax=as.name(max_error)), width=1)
  p1<-p1 + labs(x=xaxis_text, y= yaxis_text, colour="Tools") + scale_x_reverse() +scale_colour_manual(values=cbbPalette) + ggtitle(title)
  if (legend_true == FALSE){
    p1<-p1 + theme(legend.position = 'none', text = element_text(size=14), plot.title = element_text(size=14) )
  }

  return(p1)


}
###########################################

#Read in data
figure_5a_data<-read.table("../data/Figure5a.txt")

#Process data for ggplot
ES_figure_5a_data<-process_data(figure_5a_data, "ES")

#Set up colour blind friendly palette
tools<-c("RSEM","Salmon Alignment", "Salmon Quasi", "Salmon SMEM", "Sailfish", "eXpress", "Kallisto")
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(cbbPalette)<-tools

#Create ES cell graphs
ES_drop_spear<-plot_results("\nES cells", "","", FALSE, ES_figure_5a_data, "spearmans_results", "max_spearmans_error", "min_spearmans_error")
ES_drop_NRMSE<-plot_results("", "","Threshold % Dropouts", FALSE, ES_figure_5a_data, "nrmse_results", "max_nrmse_error", "min_nrmse_error")

#Create figure legends
leg<-g_legend(plot_results("Spearman's Rho", "Spearman's Rho","Threshold Percentage of Dropouts", TRUE, ES_figure_5a_data, "spearmans_results", "max_spearmans_error", "min_spearmans_error"))

#Arrange graphs
Figure5a<-ggarrange(ggarrange(ES_drop_spear + ylim(0.6,1), ES_drop_NRMSE + ylim(0,50),nrow=2, ncol=1), leg,
                    nrow = 1,
                    ncol = 2
)

rm(list=setdiff(ls(), "Figure5a"))

#######################################################
#Functions

#Function to find point density
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#Function to plot graphs
plot_cor_expr_zeros<-function(df, ID, x_lab_text, y_lab_text, title, legend){

  #filter by ID
  df<-df[df$ID==as.name(ID),]

  #find density
  df$Density <- get_density(df$expression, df$percent_zeros)

  #create graph
  p1<-ggplot(data=df, aes(x=expression, y=percent_zeros, colour=Density)) + geom_point() + scale_color_viridis()
  p1<-p1 + xlab(x_lab_text) + ylab(y_lab_text) + ggtitle(title) + theme(text = element_text(size=14), plot.title = element_text(size=14))
  if (legend==FALSE){
    p1<- p1 + theme(legend.position = 0)
  }

  return(p1)
}

# Function to create legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

#########################################################################

Figure_5b_data<-read.table("../data/Figure5b.txt")

leg_cor<-g_legend(plot_cor_expr_zeros(Figure_5b_data, "simulated_ES", " ", "% zeros", "Simulated ES cells", TRUE))

ES_sim_cor_expr_zeros<-plot_cor_expr_zeros(Figure_5b_data, "simulated_ES"," ", " ", "\nSimulated ES cells", FALSE)
ES_real_cor_expr_zeros<-plot_cor_expr_zeros(Figure_5b_data, "real_ES", "log2(counts + 1)", " ", "\nReal ES cells", FALSE)

Figure5b<-ggarrange(ggarrange(ES_sim_cor_expr_zeros, ES_real_cor_expr_zeros, ncol=1, nrow=2), leg_cor, ncol = 2, widths=c(2,1))

ggarrange(Figure5a, Figure5b, ncol=1, nrow=2)

ggsave("../pdfs/Figure5.pdf", plot=last_plot(), height= 225, width=170, units=c("mm"))
ggsave("../pngs/Figure5.png", last_plot(), height= 225, width=170, units=c("mm"))
