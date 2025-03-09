library(edgeR)
library(dplyr)
library(showtext)

table <- read.table("CFP1.RNAseq.counts.txt",header = T)
rownames(table) <- table$GeneSymbol
table <- table[,-1]
cpm_log <- cpm(table, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
expr_cutoff <- -1
table2 <- table[median_log2_cpm > expr_cutoff, ]
group <- c("WT","WT","KO","KO")
y <- DGEList(counts = table2, group = group)
y <- calcNormFactors(y,method = "TMM")
y <- estimateCommonDisp(y)
et <- exactTest(y,pair = c("WT","KO"))
results_edgeR <- topTags(et, n = nrow(table2), sort.by = "none")
result_table <- results_edgeR$table
result_table$GeneSymbol <- rownames(result_table)
counts.per.m <- cpm(y, normalized.lib.sizes=TRUE)
counts.per.m <- data.frame(counts.per.m)
counts.per.m$GeneSymbol <- rownames(table2)
counts.per.m2 <- merge(counts.per.m,result_table,by="GeneSymbol")

################################ FigS6F
up <- counts.per.m2[counts.per.m2$logFC>0.5849625 & counts.per.m2$PValue< 0.05,]
down <- counts.per.m2[counts.per.m2$logFC< -0.5849625 & counts.per.m2$PValue< 0.05,]
nonregulated <- counts.per.m2[!counts.per.m2$GeneSymbol %in% up$GeneSymbol & !counts.per.m2$GeneSymbol %in% down$GeneSymbol,]
up$type <- "Up"
down$type <- "Down"
nonregulated$type <- "Nonregulated"
allgenes <- rbind(up,down,nonregulated)
library(ggplot2)
allgenes$type <- factor(allgenes$type,levels = c("Up","Down","Nonregulated"))
ggplot(allgenes,aes(x=logFC,y=-log10(PValue),col=type)) + geom_point() + coord_cartesian(xlim = c(-6,6),ylim = c(0,60))  + theme_classic() + theme(legend.position = "none") +
  scale_color_manual(values=c("red","blue","grey"))

################################ FigS6H

library(pheatmap)
rownames(allgenes) <- allgenes$GeneSymbol
genes <- data.frame(GeneSymbol=c("Il10","Tigit","Lag3","Ccr4","Ccr6","Icos","Klrg1","Tcf7","Pdcd1","Ctla4","Nt5e","Cxcr3","Tnfrsf18","Cxcr6","Fgl2","Itgae","Il1rl1","Folr4"))
selectedgenes <- genes %>% left_join(allgenes,by="GeneSymbol")
rownames(selectedgenes) <- selectedgenes$GeneSymbol
selectedgenes <- selectedgenes[,-1]
pheatmap((selectedgenes[,c(1:4)]),scale = "row",cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("blue","white","red"))(101),border_color = NA,fontsize = 12,show_rownames = T)

################################ FigS6G

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

go_down = run_GO_mice(candidate_gene = down$GeneSymbol, background_gene = NULL)
