# Fig7C-D
# calculate enrichment levels of H3K4me3 in WT and KO cells around Broad TSS using Source Code 1

broadTSS <- read.table("H3K4me3-WTKO-BroadTSS1.overlap.sort.gene.txt")
broadTSS <- broadTSS[broadTSS$V14==0,]
broadTSS$V13 <- as.numeric(broadTSS$V13)
broadTSS$V9 <- as.numeric(broadTSS$V9)

broadTSS$FC <- log2(broadTSS$V13 + 0.1) - log2(broadTSS$V9 + 0.1)

broadTSS <- broadTSS[order(broadTSS$FC),]
broadTSS_down <- broadTSS[broadTSS$FC< -1,]
write.table(broadTSS_down,"H3K4me3-WTKO-BroadTSS1.overlap.sort.gene.down.txt",col.names=T,row.names=F,sep="\t",quote=F)


x1 <- read.table("Foxp3cre-CFP1-2reps-20230713-IgG.sort.merge.TSSpromoter2kb.bed")
x2 <- read.table("GSE121279_peaks-foxp3.mm9.TSSpromoter2kb.bed")
x3 <- read.table("H3K4me3-WTKO-BroadTSS1.overlap.sort.gene.down.txt",header = T)
library(VennDiagram)
T1<-venn.diagram(x =list(Cfp1_Promoter_Gene=(x1$V5),Foxp3_Promoter_Gene=(x2$V5),H3K4me3KOdown_BroadTSS=x3$genename), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("blue","red","yellow"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0, 0,180))
grid.draw(T1)

s <- intersect(x1$V5,x2$V5)
s <- intersect(s,x3$genename)

ego = run_GO_mice(candidate_gene = s, background_gene = NULL)
