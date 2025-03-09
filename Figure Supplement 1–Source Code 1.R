
## FigS1A
##########

x <- read.table("Treg-H3K4me3_WT_peaks.bed")
x <- x[,c(1,2,3)]
colnames(x) <- c("chr","start","end")

library(methylKit)
x2 <- as(x,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)
diffAnn = annotateWithGeneParts(as(x2,"GRanges"),gene.obj)

plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="Treg-WT-H3K4me3",cex.legend=0.5)


########## obtain random regions using bedtools

x <- read.table("Treg-H3K4me3_WT_peaks.shuffle.bed")
x <- x[,c(1,2,3)]
colnames(x) <- c("chr","start","end")

library(methylKit)
x2 <- as(x,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)
diffAnn = annotateWithGeneParts(as(x2,"GRanges"),gene.obj)

plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="random H3K4me3",cex.legend=0.5)

