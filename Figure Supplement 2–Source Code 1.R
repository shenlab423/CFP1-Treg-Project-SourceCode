####################### FigS2B

x <- read.table("Foxp3cre-CFP1-2reps-20230713-IgG.sort.merge.bed")
x <- x[,c(1,2,3)]
colnames(x) <- c("chr","start","end")

library(methylKit)
x <- as(x,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)

diffAnn2 = annotateWithGeneParts(as(x,"GRanges"),gene.obj)


plotTargetAnnotation(diffAnn2,precedence=TRUE,
                     main="CFP1 peaks annotation",cex.legend=0.5)



x <- read.table("Foxp3cre-CFP1-2reps-20230713-IgG.sort.merge.shuffle.bed")
x <- x[,c(1,2,3)]
colnames(x) <- c("chr","start","end")

library(methylKit)
x <- as(x,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)

diffAnn2 = annotateWithGeneParts(as(x,"GRanges"),gene.obj)


plotTargetAnnotation(diffAnn2,precedence=TRUE,
                     main="Random peaks annotation",cex.legend=0.5)

####################### FigS2D

library(ChIPpeakAnno)

x1 <- read.table("Foxp3cre-CFP1-2reps-20230713-IgG.sort.merge.bed")
x1 <- x1[,c(1,2,3)]
colnames(x1) <- c("seqnames","start","end")
x1 <- unique(x1)
Cfp1 <- toGRanges(x1)


x2 <- read.table("GSE121279_peaks-foxp3.mm9.bed")
x2 <- x2[,c(1,2,3)]
colnames(x2) <- c("seqnames","start","end")
x2 <- unique(x2)
Foxp3 <- toGRanges(x2)


x <- findOverlapsOfPeaks(Cfp1,Foxp3,connectedPeaks="keepAll")
makeVennDiagram(x,connectedPeaks="keepAll")