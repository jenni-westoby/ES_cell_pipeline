library(tidyverse)

###########################################################
# FUNCTIONS

make_df<-function(df, name_df, origin_df){
  percentage_df<-sum(df$tpm==0)/nrow(df) * 100
  bulk_df<-data.frame(cell=name_df, percentage=percentage_df, origin=origin_df)
  return(bulk_df)
}

make_zeros_plot<-function(df){
  graph<-ggplot(data=df, aes(y=percentage, x=origin)) + geom_jitter(alpha=0.5, position=position_jitter(width = .2))
  graph<-graph + ylab("% of Isoforms Which Are Unexpressed") 
}

########################


#read in Kolod et al sc #zeros data
sc_kolod<-read.table(gzfile("../data/SupplementaryFigure12.gz"))


sc_kolod<-sc_kolod %>% group_by(cell) %>% summarise((sum(estimates==0)/128809)*100)
names(sc_kolod)[2]<-c("percentage")
sc_kolod<-data.frame(sc_kolod, origin="ES single cell")


#bulk Kolod et al expression matrix
bulk_ES<-read.table("../../Simulation/Kallisto_results_real_data/ERR522956/abundance.tsv", header=T)
bulk_ES<-make_df(bulk_ES, "bulk_ES", "ES bulk")

master_df<-rbind(sc_kolod, bulk_ES)
make_zeros_plot(master_df)

ggsave("../pdfs/SupplementaryFigure17.pdf", plot=last_plot())
ggsave("../pngs/SupplementaryFigure17.png", plot=last_plot())
