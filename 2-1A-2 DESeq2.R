## Install packages
install.packages(c("gplots","genefilter","calibrate","RColorBrewer","DESeq2"))

## Load required libraries
library("DESeq2")
library("RColorBrewer")
library("calibrate")
library("genefilter")
library("gplots")

## Combine count files into dataframe
# Import data from featureCounts
countdata <- read.table("counts.txt", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove sorted.bam from filenames
colnames(countdata) <- gsub("\\sorted.bam$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition 
sampleCondition <- c("SI","Colon","Lung","MWAT","SI","Colon","Lung","MWAT","SI","Colon","Lung","MWAT")

# Analysis with DESeq2 ----------------------------------------------------
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), sampleCondition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~sampleCondition)

## Run DESeq normalization
dds<-DESeq(ddsHTSeq)

## Plot dispersions
svg("Dispersion Plot.svg", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion Plot")
dev.off()

## Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

## Sample distance heatmap
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(sampleCondition))])
sampleDists <- as.matrix(dist(t(assay(rld))))
svg("Sample Distance Heatmap.svg", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[sampleCondition], RowSideColors=mycols[sampleCondition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

## Principal components analysis
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) 
  {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
svg("PCA.svg", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()

## Differential Expression (multiple tissues)
res.SI.Colon <- results(dds, contrast=c("condition","SI","Colon"))
res.SI.Lung <- results(dds, contrast=c("condition","SI","Lung"))
res.SI.MWAT <- results(dds, contrast=c("condition","SI","MWAT"))
res.Colon.Lung <- results(dds, contrast=c("condition","Colon","Lung"))
res.Colon.MWAT <- results(dds, contrast=c("condition","Colon","MWAT"))
res.Lung.MWAT <- results(dds, contrast=c("condition","Lung","MWAT"))

## Filter results with padj < 0.05
#table(res.SI.Colon$padj<0.5)
#table(res.SI.Lung$padj<0.5)
#table(res.SI.MWAT$padj<0.5)
#table(res.Colon.Lung$padj<0.5)
#table(res.Colon.MWAT$padj<0.5)
#table(res.Lung.MWAT$padj<0.5)

#Order by P-value
res.SI.Colon <- res.SI.Colon[order(res.SI.Colon$padj),]
res.SI.Lung <- res.SI.Lung[order(res.SI.Lung$padj),]
res.SI.MWAT <- res.SI.MWAT[order(res.SI.MWAT$padj),]
res.Colon.Lung <- res.Colon.Lung[order(res.Colon.Lung$padj),]
res.Colon.MWAT <- res.Colon.MWAT[order(res.Colon.MWAT$padj),]
res.Lung.MWAT <- res.Lung.MWAT[order(res.Lung.MWAT$padj),]

## Write DE results
write.csv(res.SI.Colon,file="res.SI.Colon-DE.csv")
write.csv(res.SI.Lung,file="res.SI.Lung-DE.csv")
write.csv(res.SI.MWAT,file="res.SI.MWAT-DE.csv")
write.csv(res.Colon.Lung,file="res.Colon.Lung-DE.csv")
write.csv(res.Colon.MWAT,file="res.Colon.MWAT-DE.csv")
write.csv(res.Lung.MWAT,file="res.Lung.MWAT-DE.csv")

## Examine plot of p-values
##hist(res$pvalue, breaks=50, col="grey")

## MA plot
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) 
  {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
svg("res.SI.Colon-MA.svg", 1500, 1000, pointsize=20)
maplot(res.SI.Colon, main="MA Plot")
dev.off()
svg("res.SI.Lung-MA.svg", 1500, 1000, pointsize=20)
maplot(res.SI.Lung, main="MA Plot")
dev.off()
svg("res.SI.MWAT-MA.svg", 1500, 1000, pointsize=20)
maplot(res.SI.MWAT, main="MA Plot")
dev.off()
svg("res.Colon.Lung-MA.svg", 1500, 1000, pointsize=20)
maplot(res.Colon.Lung, main="MA Plot")
dev.off()
svg("res.Colon.MWAT-MA.svg", 1500, 1000, pointsize=20)
maplot(res.Colon.MWAT, main="MA Plot")
dev.off()
svg("res.Lung.MWAT-MA.svg", 1500, 1000, pointsize=20)
maplot(res.Lung.MWAT, main="MA Plot")
dev.off()



## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) 
  {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

svg("res.SI.Colon-VP.svg", 1200, 1000, pointsize=20)
volcanoplot(res.SI.Colon, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), main="Volcano Plot")
dev.off()
svg("res.SI.Lung-VP.svg", 1200, 1000, pointsize=20)
volcanoplot(res.SI.Lung, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), main="Volcano Plot")
dev.off()
svg("res.SI.MWAT-VP.svg", 1200, 1000, pointsize=20)
volcanoplot(res.SI.MWAT, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), main="Volcano Plot")
dev.off()
svg("res.Colon.Lung-VP.svg", 1200, 1000, pointsize=20)
volcanoplot(res.Colon.Lung, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), main="Volcano Plot")
dev.off()
svg("res.Colon.MWAT-VP.svg", 1200, 1000, pointsize=20)
volcanoplot(res.Colon.MWAT, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), main="Volcano Plot")
dev.off()
svg("res.Lung.MWAT-VP.svg", 1200, 1000, pointsize=20)
volcanoplot(res.Lung.MWAT, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2), main="Volcano Plot")
dev.off()



##topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
##heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" ), [colData(rld)$condition ] )


