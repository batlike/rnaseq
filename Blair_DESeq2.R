# Initiate Bioconductor
source("http://bioconductor.org/biocLite.R")

# Load Libraries and Import Data
library('DESeq2')

# Import data from htseq-counts
directory<-"Y:/RUShares/Lab Members/Blair/Harvard RNASeq Data/Counts/test"
sampleFiles <- grep("ctrl",list.files(directory),value=TRUE)

# Assign sample condition to dataframe
sampleCondition<-c("ctrl rt","ctrl rt","ctrl rt","ctrl cold","ctrl cold","ctrl cold")
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

# Load DESeq dataset from HTSeqCount 
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)

# Load conditions
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("ctrl rt","ctrl cold"))

# DESeq normalization
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)

# MA-Plot
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.copy(png,"deseq2_MAplot.png")
dev.off()
###


