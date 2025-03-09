######## Fig1A-B


library(ChIPpeakAnno)

x1 <- read.table("Treg-H3K4me3_WT_peaks.bed")
x1 <- x1[,c(1,2,3)]
colnames(x1) <- c("seqnames","start","end")
x1 <- unique(x1)
H3K4me3 <- toGRanges(x1)


x2 <- read.table("GSE121279_peaks-foxp3.mm9.bed")
x2 <- x2[,c(1,2,3)]
colnames(x2) <- c("seqnames","start","end")
x2 <- unique(x2)
Foxp3 <- toGRanges(x2)


x <- findOverlapsOfPeaks(H3K4me3,Foxp3,connectedPeaks="keepAll")
makeVennDiagram(x,connectedPeaks="keepAll")


H3K4me3_foxp3 <- data.frame(x$peaklist$`H3K4me3///Foxp3`)

H3K4me3_foxp3 <- H3K4me3_foxp3[,c(1,2,3)]
colnames(H3K4me3_foxp3) <- c("chr","start","end")

write.table(H3K4me3_foxp3,"H3K4me3_foxp3_overlapped.bed",col.names = F,row.names = F,sep="\t",quote=F)

library(methylKit)
H3K4me3_foxp3 <- as(H3K4me3_foxp3,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)
H3K4me3_foxp3 = annotateWithGeneParts(as(H3K4me3_foxp3,"GRanges"),gene.obj)

plotTargetAnnotation(H3K4me3_foxp3,precedence=TRUE,
                     main="H3K4me3 Foxp3 overlapped peaks",cex.legend=0.5)


########


library(ChIPpeakAnno)

x1 <- read.table("nTreg-H3K27me3.q5_peaks.broadPeak")
x1 <- x1[,c(1,2,3)]
colnames(x1) <- c("seqnames","start","end")
x1 <- unique(x1)
H3K27me3 <- toGRanges(x1)


x2 <- read.table("GSE121279_peaks-foxp3.mm9.bed")
x2 <- x2[,c(1,2,3)]
colnames(x2) <- c("seqnames","start","end")
x2 <- unique(x2)
Foxp3 <- toGRanges(x2)


x <- findOverlapsOfPeaks(H3K27me3,Foxp3,connectedPeaks="keepAll")
makeVennDiagram(x)

H3K27me3_foxp3 <- data.frame(x$peaklist$`H3K27me3///Foxp3`)

H3K27me3_foxp3 <- H3K27me3_foxp3[,c(1,2,3)]
colnames(H3K27me3_foxp3) <- c("chr","start","end")

library(methylKit)
H3K27me3_foxp3 <- as(H3K27me3_foxp3,"GRanges")

library(genomation)
gene.obj=readTranscriptFeatures("mm9.bed",up.flank=2000,down.flank=2000)
H3K27me3_foxp3 = annotateWithGeneParts(as(H3K27me3_foxp3,"GRanges"),gene.obj)

plotTargetAnnotation(H3K27me3_foxp3,precedence=TRUE,
                     main="H3K27me3 Foxp3 overlapped peaks",cex.legend=0.5)



