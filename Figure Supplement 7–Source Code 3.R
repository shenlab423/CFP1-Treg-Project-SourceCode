########## FigS7E

library(ChIPpeakAnno)

x1 <- read.table("Treg-H3K4me1_WT_peaks.bed")
x1 <- x1[,c(1,2,3)]
colnames(x1) <- c("seqnames","start","end")
x1 <- unique(x1)
WTH3K4me1 <- toGRanges(x1)


x2 <- read.table("Treg-H3K4me1_KO_peaks.bed")
x2 <- x2[,c(1,2,3)]
colnames(x2) <- c("seqnames","start","end")
x2 <- unique(x2)
KOH3K4me1 <- toGRanges(x2)


x <- findOverlapsOfPeaks(WTH3K4me1,KOH3K4me1,connectedPeaks="keepAll")
makeVennDiagram(x)


########## FigS7F

x <- read.table("Treg-H3K4me1_WT_peaks.bed")
x <- x[,c(1,2,3)]
colnames(x) <- c("chr","start","end")

library(methylKit)
x2 <- as(x,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)
diffAnn = annotateWithGeneParts(as(x2,"GRanges"),gene.obj)

plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="Treg-WT-H3K4me1",cex.legend=0.5)


##########

x <- read.table("Treg-H3K4me1_KO_peaks.bed")
x <- x[,c(1,2,3)]
colnames(x) <- c("chr","start","end")

library(methylKit)
x2 <- as(x,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)
diffAnn = annotateWithGeneParts(as(x2,"GRanges"),gene.obj)

plotTargetAnnotation(diffAnn,precedence=TRUE,
                     main="Treg-KO-H3K4me1",cex.legend=0.5)

