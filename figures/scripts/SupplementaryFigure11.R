#QC script
library(ggplot2)
library(ggpubr)
library(scater, quietly = TRUE)
library(knitr)
options(stringsAsFactors = FALSE)

############################################################
# READS AND ALIGNMENT QC PLOTS

QC_raw<-read.table("../data/SupplementaryFigure11_reads_alignment_data.txt")

#Alignment based QC plots
Unique<-ggplot(data=QC_raw, aes(x=reorder(Filename,Unique), y=Unique/1000000)) + geom_point(stat="identity") + ylab("Millions of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() +xlab(" ") + ggtitle("\nNumber of Uniquely Mapping Reads")
Unique<-Unique + geom_hline(yintercept=5000000/1000000, color='red', linetype='dashed')

NonUnique<-ggplot(data=QC_raw, aes(x=reorder(Filename,NonUnique), y=NonUnique/1000000)) + geom_point(stat="identity") + ylab("Millions of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + xlab(" ") + ggtitle("Number of Non-Uniquely Mapping Reads")
#NonUnique<-NonUnique + geom_hline(yintercept=2000000/1000000, color='red', linetype='dashed')

Unmapped<-ggplot(data=QC_raw, aes(x=reorder(Filename,Unmapped), y=Unmapped)) + geom_point(stat="identity") + ylab("Number of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + ylim(0,100) + xlab(" ") + ggtitle("Number of Unmapped Reads")
NumAlign<-ggplot(data=QC_raw, aes(x=reorder(Filename,NumAlignments), y=NumAlignments/1000000)) + geom_point(stat="identity") + ylab("Millions of Alignments") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_x_discrete() + xlab(" ") + ggtitle("\nNumber of Alignments")

#Read number based QC plot
NumReads<-ggplot(data=QC_raw, aes(x=reorder(Filename,NumReads), y=NumReads/1000000)) + geom_point(stat="identity") + ylab("Millions of Reads") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_text(size=14)) + scale_x_discrete() + xlab(" ") + ggtitle("\nNumber of Reads")

#############################################################
# MT READS QC PLOT

load("../data/SupplementaryFigure11_scater_object.RData")

mt_reads<-plotPhenoData(
  teich_scater_QC,
  aes_string(x = "total_features",
             y = "pct_counts_MT",
             colour = "batch")
)

# make final figure
ggarrange(mt_reads, NumReads, NumAlign, Unique, NonUnique, Unmapped, ncol = 2,nrow=3, labels = c("A","B","C","D","E","F"))
ggsave("../pdfs/SupplementaryFigure11.pdf", plot=last_plot(),height = (6.04*1.5), width=(8.57*1.5))
ggsave("../pngs/SupplementaryFigure11.png", plot=last_plot(),height = (6.04*1.5), width=(8.57*1.5))