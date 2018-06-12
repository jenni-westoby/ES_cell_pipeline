library(ggplot2)
library(ggpubr)

#read in data
ground_truth_expr<-read.table(gzfile("../data/SupplementaryFigure12.gz"))

#plot data and save
ggplot(data=ground_truth_expr, aes(x=estimates)) + geom_histogram(binwidth = 10)  + coord_cartesian(ylim = c(0, 100000)) + xlab("100 x log2(Ground Truth Expression in TPM + 1)") + ylab("Frequency of Isoforms") 
ggsave("../pdfs/SupplementaryFigure12.pdf", plot=last_plot())
ggsave("../pngs/SupplementaryFigure12.png", plot=last_plot())