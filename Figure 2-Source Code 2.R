#### Histone modification enrichment was calculated using Source Code 1
#### Fig2B 

promoter <- read.table("mm9_promoter.bed")

xtable1 <- read.table("H3K4me1/H3K4me1-WT-0520-mm9_promoter.rpkm")
xtable3 <- read.table("H3K4me1/07221-H3K4me1-WT-mm9_promoter.rpkm")

xtable5 <- read.table("H3K4me3/H3K4-WT-mm9_promoter.rpkm")
xtable7 <- read.table("H3K4me3/H3K4me3-WT-0520-mm9_promoter.rpkm")

xtable <- cbind(promoter,xtable1$V4,xtable3$V4,xtable5$V4,xtable7$V4)
colnames(xtable) <- c("chr","start","end","NM","GeneName","H3K4me1-WT1","H3K4me1-WT2","H3K4me3-WT1","H3K4me3-WT2")

library(dplyr)
xtable <- xtable %>% group_by(GeneName) %>% summarise(`H3K4me1-WT1`=mean(`H3K4me1-WT1`),`H3K4me1-WT2`=mean(`H3K4me1-WT2`),
                                                      `H3K4me3-WT1`=mean(`H3K4me3-WT1`),`H3K4me3-WT2`=mean(`H3K4me3-WT2`))

xtable2 <- read.table("nTreg-H3K27me3-mm9_promoter.rpkm")
colnames(xtable2) <- c("chr","start","end","NM","GeneName","nTreg-H3K27me3")
xtable2 <- xtable2 %>% group_by(GeneName) %>% summarise(`nTreg-H3K27me3`=mean(`nTreg-H3K27me3`))


xtable <- cbind(xtable,xtable2)
all <- xtable[,c(1,2,3,4,5,7)]
all$`H3K4me1-WT` <- 0.5 * (all$`H3K4me1-WT1` + all$`H3K4me1-WT2`)
all$`H3K4me3-WT` <- 0.5 * (all$`H3K4me3-WT1` + all$`H3K4me3-WT2`)

K43enrich <- all[all$`H3K4me3-WT`>8,]

## obtain FOXP3 or CFP1 peaks covered genes using bedtools (from TSS-2kb to TES)

foxp3 <- read.table("GSE121279_peaks-foxp3.mm9.TSSpromoter2kb.bed")
cfp1 <- read.table("Foxp3cre-CFP1-2reps-20230713-IgG.sort.merge.TSSpromoter2kb.bed")

library(VennDiagram)

T1<-venn.diagram(x =list(Foxp3=foxp3$V5,Cfp1=cfp1$V5,H3K4me3=K43enrich$GeneName), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("blue","red","yellow"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0, 0, 180))
grid.draw(T1)

