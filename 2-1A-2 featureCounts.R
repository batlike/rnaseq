biocLite("Rsamtools")
biocLite("Rsubread")
library(Rsubread)
library(Rsamtools)

## Sort BAM files and generate index .bai files

types <- c("SI", "Colon", "Lung", "MWAT")
samples <- c("1", "2", "3")

sortIndexBam <- function (types,samples){
for (type in types)
  {
  for (sample in samples)
  {
    sortBam(
      paste0(type, sample, ".bam"),
      paste0(type, sample, "sorted"))
    indexBam(
      paste0(type, sample, "sorted"))
  }
  }
}
sortIndexBam(types,samples)


## Get and export feature counts from aligned BAM files

counts <- featureCounts("SI1sorted.bam",
                        "Colon1sorted.bam",
                        "Lung1sorted.bam",
                        "MWAT1sorted.bam",
                        "SI2sorted.bam",
                        "Colon2sorted.bam",
                        "Lung2sorted.bam",
                        "MWAT2sorted.bam",
                        "SI3sorted.bam",
                        "Colon3sorted.bam",
                        "Lung3sorted.bam",
                        "MWAT3sorted.bam",
              annot.ext="Mus_musculus.GRCm38.81.gtf",
              isGTFAnnotationFile=TRUE,
              GTF.featureType="exon",
              GTF.attrType="gene_id",
              useMetaFeatures=TRUE,
              allowMultiOverlap=FALSE,
              isPairedEnd=FALSE,
              requireBothEndsMapped=FALSE,
              checkFragLength=FALSE,
              minFragLength=50,
              maxFragLength=600,
              nthreads=1,
              strandSpecific=0,
              minMQS=0,
              readExtension5=0,
              readExtension3=0,
              read2pos=NULL,
              minReadOverlap=1,
              countSplitAlignmentsOnly=FALSE,
              countMultiMappingReads=FALSE,
              countPrimaryAlignmentsOnly=FALSE,
              countChimericFragments=TRUE,
              ignoreDup=FALSE,
              chrAliases=NULL,
              reportReads=FALSE)

write.table(x=data.frame(counts$annotation[,c("GeneID")],
                         counts$counts,
                         stringsAsFactors=FALSE),
            file="counts.txt",
            quote=FALSE,sep="\t",
            row.names=FALSE)

