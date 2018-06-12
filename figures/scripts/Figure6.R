library(ggplot2)
library(ggpubr)
library(reshape2)

##############################
# FUNCTIONS

get_percentage_hist<-function(df){
  no_zeros<-df#[df$value!=0,]
  no_zeros$value<-no_zeros$value==2
  two_iso_percentage<-no_zeros %>% group_by(genes) %>% summarise(percentage(value))
  ggplot(data = two_iso_percentage, aes(x=`percentage(value)`)) + geom_histogram() + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")
}

#####################################
# MAKE FIGURE 6A

barplot_data<-read.table("../data/Figure6_number_of_isoforms.txt")
ES_barplot_data<-barplot_data[barplot_data$ID=="ES",]

ES_0_1_2_barplot<-ggplot(data=ES_barplot_data, aes(x=value)) + geom_bar() + xlab("Number of Isoforms Expressed") + ggtitle("ES cells")

##############################################################
#MAKE FIGURE 6B

exprs_percent_data<-read.table("../data/Figure_6_percent_exprs.txt")
colnames(exprs_percent_data)[2]<-"expression"
colnames(exprs_percent_data)[3]<-"percentage"

ES_exprs_percent_data<-exprs_percent_data[exprs_percent_data$ID=="ES",]

ES_percent_hist<-ggplot(data=ES_exprs_percent_data, aes(x=percentage)) + geom_histogram() + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")

#################################################################
#MAKE FIGURE 6C

ES_exprs_percent_plot<-ggplot(data = ES_exprs_percent_data, aes(y=expression,x=percentage)) + geom_point() + ylab("log2(Counts + 1)") + xlab("% Cells Which Express Both Isoforms") + ggtitle(" ")

ggarrange(ES_0_1_2_barplot, ES_percent_hist, ES_exprs_percent_plot, ncol =1, nrow=3, labels=c("A", "B", "C"))
ggsave("../pdfs/Figure6.pdf", plot=last_plot(), height= 225, width=170, units=c("mm"))
ggsave("../pngs/Figure6.png", plot=last_plot(), height= 225, width=170, units=c("mm"))
