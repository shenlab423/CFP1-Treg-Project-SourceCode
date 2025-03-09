# Fig1D-F
#### Histone modification enrichment was calculated using Source Code 1

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


all <- all[,c(1,8,6)]
set.seed(1000)
km <- kmeans(all[,2:3],5)
km

all4 <- cbind(all,km$cluster)

all4 <- all4[order(all4$`km$cluster`),]

all4$`H3K4me3-WT` <- log2(all4$`H3K4me3-WT`+0.1)
all4$`nTreg-H3K27me3` <- log2(all4$`nTreg-H3K27me3`+0.1)

mm9 <- read.table("mm9_promoter.bed")
colnames(mm9) <- c("chr","start","end","NM","GeneName")
library(dplyr)
x4 <- all4 %>% left_join(mm9,by="GeneName")
x4 <- x4[,c(5,6,7,8,1)]
write.table(x4,"HM_clustered_promoter.nTreg.bed",col.names = F,row.names = F,sep = "\t",quote=F)

# Heatmap was generated using Figure 1-Source Code 3

###########################
x <- read.table("genes.ref.fpkm_table",header = T)
orders <- read.table("HM_clustered_promoter.nTreg.bed")

library(dplyr)

colnames(x) <- c("V5","WT1","WT2","WT3","WT4")

x$WTave <- 0.25 * (x$WT1 + x$WT2 + x$WT3 + x$WT4)

x2 <- orders %>% left_join(x,by="V5")

x2 <- x2[,c(5,5,6,7,8,9,10)]

x2$WT1 <- log2(x2$WT1 + 0.1)
x2$WT2 <- log2(x2$WT2 + 0.1)
x2$WT3 <- log2(x2$WT3 + 0.1)
x2$WT4 <- log2(x2$WT4 + 0.1)
x2$WTave <- log2(x2$WTave + 0.1)

write.table(x2,"nTreg.ref.RNA-seq.txt",col.names = T,row.names = F,sep="\t",quote=F)

## visualize Fig1E using JavaTreeView with file "nTreg.ref.RNA-seq.txt"
###########################

library(ggplot2)
library(ggpubr)
library(tibble)
library(data.table)
library(DT)
library(org.Mm.eg.db)
library(clusterProfiler, quietly = T)
library(DOSE, quietly = T)

set.seed(1000)

run_GO_mice = function(candidate_gene, background_gene=NULL, gene_format = 'SYMBOL', ontology = 'BP', cutoff=0.05,
                       showCategory=10,font.size=10,title = 'GO enrichment'){
  diff_gene_ID<-clusterProfiler::bitr(candidate_gene, fromType = gene_format, toType ="ENTREZID", OrgDb="org.Mm.eg.db")
  if(is.null(background_gene)){
    ego <-  simplify(enrichGO(gene = diff_gene_ID$ENTREZID,  OrgDb = org.Mm.eg.db,
                              keyType = 'ENTREZID', ont = ontology, readable = T,
                              pAdjustMethod = "BH", qvalueCutoff  = cutoff, pvalueCutoff  = cutoff))
  } else{
    background_gene = clusterProfiler::bitr(background_gene, fromType = gene_format, toType ="ENTREZID", OrgDb="org.Mm.eg.db")
    ego <-  simplify(enrichGO(gene = diff_gene_ID$ENTREZID,  OrgDb = org.Mm.eg.db,
                              universe = background_gene$ENTREZID,
                              keyType = 'ENTREZID', ont = ontology, readable = T,
                              pAdjustMethod = "BH", qvalueCutoff  = cutoff, pvalueCutoff  = cutoff))
  }
  
  if(nrow(ego@result)>0){
    print(dotplot(ego, showCategory = showCategory, font.size = font.size, x='Count',title=title))
  }
  return(ego)
}

ego = run_GO_mice(candidate_gene = all4$GeneName[all4$`km$cluster`<3], background_gene = NULL)
ego = run_GO_mice(candidate_gene = all4$GeneName[all4$`km$cluster`==3], background_gene = NULL)
ego = run_GO_mice(candidate_gene = all4$GeneName[all4$`km$cluster`==4], background_gene = NULL)
ego = run_GO_mice(candidate_gene = all4$GeneName[all4$`km$cluster`==5], background_gene = NULL)

###########################
