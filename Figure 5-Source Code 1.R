### Note 
# The FASTQ files were processed using the BD Rhapsody Targeted Analysis Pipeline
# on the Seven Bridges platform (https://www.sevenbridges.com). Following processing, 
# we obtained the Foxp3cre-TCR_1_Seurat.rds file, which was used for downstream analysis.

# V(D)J reads were processed with the BD Rhapsody V(D)J CDR3 Pipeline on Seven Bridges, 
# resulting in the Foxp3cre-TCR_1_VDJ_perCell.csv file, which was also used for downstream analysis.

###


### Load required packages

library(Seurat)
library(harmony)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(slingshot)
library(SingleCellExperiment)
library(colorRamps)
library(pheatmap)
library(scRepertoire)


tregs = readRDS("Foxp3cre-TCR_1_Seurat.rds")
tregs = subset(tregs, Sample_Tag %in% c("Multiplet", "Undetermined"), invert = T)

tregs <- NormalizeData(tregs)
all.genes <- rownames(tregs)
tregs <- ScaleData(tregs, vars.to.regress = "nCount_RNA", features = all.genes)

tregs <- RunPCA(tregs, features = all.genes, verbose = FALSE)
ElbowPlot(tregs, ndims = 50)
p1 <- DimPlot(tregs, reduction = "pca", group.by = "Sample_Tag")

tregs <- RunHarmony(tregs, "Sample_Tag", plot_convergence = TRUE)
p2 <- DimPlot(tregs, reduction = "harmony", group.by = "Sample_Tag")

cowplot::plot_grid(p1,p2)

tregs <- tregs %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1) %>% 
  identity()

DimPlot(tregs, label = T)


tregs$celltype = "NA"
tregs$celltype[tregs$seurat_clusters %in% c(4,6,10)] = "Naive subsets"
tregs$celltype[tregs$seurat_clusters == 7] = "Mki67 hi subsets"
tregs$celltype[tregs$seurat_clusters %in% c(2, 3)] = "pre actived subsets"
tregs$celltype[tregs$seurat_clusters == 9] = "Klrg1 hi actived subsets"
tregs$celltype[tregs$seurat_clusters == 1] = "Il10 hi actived subsets"
tregs$celltype[tregs$seurat_clusters == 0] = "Tnfsf14 hi actived subsets"
tregs$celltype[tregs$seurat_clusters == 5] = "Lap3 hi actived subsets"
tregs$celltype[tregs$seurat_clusters == 8] = "H2-Eb1 hi actived subsets"


index = c("Naive subsets", "pre actived subsets", 
          "Klrg1 hi actived subsets", "Il10 hi actived subsets", "Tnfsf14 hi actived subsets",
          "Lap3 hi actived subsets", "H2-Eb1 hi actived subsets", "Mki67 hi subsets")

tregs$celltype = factor(tregs$celltype, levels = index)
tregs$condition = ifelse(tregs$Sample_Tag %in% c("SampleTag01_mm", "SampleTag02_mm"), "WT", "KO")
tregs$condition = factor(tregs$condition, levels = c("WT", "KO"))


cols = c("#7a9894","#B4C5EA", "#66c2a5","#FEA891","#D6EDC3","#5AB15C",
         "#0b8F7A", "#FFD28A")

##----------------------------Figure5 A-----------------------------------------##
DimPlot(tregs, group.by = "celltype", cols = cols, pt.size = 0.2)

##---------------------------Supplementary Figure5 B----------------------------##
DimPlot(tregs, group.by = "celltype", cols = cols, pt.size = 0.2, split.by = "condition")



##----------------------------Figure5 B-----------------------------------------##

df2 <- as.data.frame.matrix(table(tregs$celltype, tregs$condition))
df2$old_per <- df2$KO/sum(df2$KO)
df2$young_per <- df2$WT/sum(df2$WT)
df2$ratio <- df2$old_per/df2$young_per
df2$log2fc <- log2(df2$ratio)


plotdat <- data.frame(celltype=rownames(df2), log2_ratio=df2$log2fc)
plotdat$celltype <- factor(plotdat$celltype, levels = index)


ggplot(plotdat, aes(x=reorder(celltype, -log2_ratio), y=log2_ratio, fill=celltype)) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  theme_bw() +  scale_fill_manual(values = cols) +
  labs(y="Log2 fold change(KO/WT)\n", x="") +
  theme_classic() + coord_flip() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    legend.title = element_blank(),
    legend.position = "NULL"
  )



##---------------------------Supplementary Figure5 C----------------------------##

df = as.data.frame(table(tregs$celltype, tregs$condition))
colnames(df)[1:2] = c("cluster", "sample")
df$cluster <- as.character(df$cluster)
df$sample = as.character(df$sample)
df$sample = factor(df$sample, levels = c("WT", "KO"))
df$cluster = factor(df$cluster, levels = index)

ggplot(df, aes(sample, Freq, fill=cluster)) +
  geom_bar(stat="identity", position="fill", width =0.7) +
  labs(y="Frequency\n", x="") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black"),
    legend.title = element_blank()
  )

##--------------------------Figure5 C------------------------------------------##

tregs_filter <- subset(tregs, celltype %in% c("Mki67 hi subsets"), invert=T)

sce <- as.SingleCellExperiment(tregs_filter)
sce <- slingshot(sce, clusterLabels = "celltype", 
                 reducedDim = "UMAP", 
                 extend = "n",
                 stretch = 0, 
                 shrink.method = "density", 
                 start.clus = "Naive subsets",
                 end.clus = c("Il10 hi actived subsets", "H2-Eb1 hi actived subsets"))

pt <- DimPlot(tregs_filter, group.by = "celltype", cols = cols, pt.size = 0.2)
pbuid <- ggplot_build(pt)
pdata <- pbuid$data[[1]]


plot(reducedDims(sce)$UMAP, col=pdata$colour, pch=16, cex=0.4)
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')



##--------------------------Supplementary Figure5 A----------------------------##

treg_heat <- read.table("treg_heatmap_gene.txt", header = T, sep = "\t")

expdat <- GetAssayData(tregs, slot = "data")
expdat <- as.matrix(expdat)
rownames(expdat) = stringr::str_split(rownames(expdat), "-", simplify = T)[,1]

expdat.split <- split(as.data.frame(t(expdat)), tregs$celltype)
plotdat2 <- Reduce(cbind, lapply(expdat.split, apply, 2, mean)) 
colnames(plotdat2) <- levels(tregs$celltype)

plotdat_checkpoint <- plotdat2[treg_heat$gene, ]

pheatmap(plotdat_checkpoint,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         cluster_cols = F,
         cluster_rows = F,
         width = 3, 
         filename = "tregs_markers_heatmap.pdf",
         border_color = "transparent",
         fontsize_row = 8,
         treeheight_row = 0)


##--------------------------Figure5 F------------------------------------------##

treg_heat <- read.table("treg_related_genes.txt", header = T, sep = "\t")
tregs <- readRDS("mxy_tregs_cluster.rds")

expdat <- GetAssayData(tregs, layer = "data")
expdat <- as.matrix(expdat)
expdat <- expdat[treg_heat$gene, ]
rownames(expdat) <- stringr::str_split(rownames(expdat), "-", simplify = T)[,1]

tregs$celltype <- paste0(tregs$condition, "_", tregs$celltype)

index = c("Naive subsets", "pre actived subsets", 
          "Klrg1 hi actived subsets", "Il10 hi actived subsets", "Tnfsf14 hi actived subsets",
          "Lap3 hi actived subsets", "H2-Eb1 hi actived subsets", "Mki67 hi subsets")
idx1 = paste0("WT_", index)
idx2 = paste0("KO_", index)

tregs$celltype <- factor(tregs$celltype, levels = c(idx1, idx2))

expdat.split <- split(as.data.frame(t(expdat)), tregs$celltype)
plotdat2 <- Reduce(cbind, lapply(expdat.split, apply, 2, mean)) 
colnames(plotdat2) <- levels(tregs$celltype)


pheatmap(plotdat2,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
         cluster_cols = F,
         #cluster_rows = T,
         width = 5, 
         filename = "tregs_markers.pdf",
         border_color = "transparent",
         fontsize_row = 8,
         treeheight_row = 0)


##--------------------------Figure5 F------------------------------------------##


dat <- processbdTCR(vdj)

contiglist <- split(dat, dat$sample_name)
samples_dat <- names(contiglist)
combined <- combineTCR(contiglist, samples = samples_dat, cells ="T-AB", filterMulti=T)

li_2 <- lapply(seq_along(combined), function(i) {
  ids <- combined[[i]][["barcode"]]
  ids1 <- sapply(strsplit(ids,"_"), "[", 3)
  combined[[i]][["barcode"]] <- ids1
  combined[[i]]})


combined2 <- do.call(rbind, li_2)
combined2$celltype <- tregs$celltype[match(combined2$barcode, colnames(tregs))]
combined2$condition <- tregs$condition[match(combined2$barcode, colnames(tregs))]

wt = combined2[combined2$condition=="WT", ]
ko = combined2[combined2$condition=="KO", ]

wt <- split(wt, wt$celltype)
ko <- split(ko, ko$celltype)

wt = clonalOverlap(wt, cloneCall="nt", method="overlap", exportTable = T)
ko = clonalOverlap(ko, cloneCall="nt", method="overlap", exportTable = T)


ko <- suppressMessages(reshape2::melt(ko))[, -1]
wt <- suppressMessages(reshape2::melt(wt))[, -1]

ko$group <- "KO"
wt $group <- "WT"

df <- rbind( wt, ko)
df$group <- factor(df$group, levels = c("WT",  "KO"))

ggplot(df, aes(x = names, y = variable, fill = value)) + 
  geom_tile() + labs(fill = "overlap") + 
  geom_text(aes(label = round(value, digits = 3))) + 
  theme_classic() + facet_wrap(~group) +
  # scale_fill_gradient(low = "white", high = 'red', na.value = "white") +
  scale_fill_gradient2(high = '#cf5246',
                       mid = 'white',
                       midpoint = ((range(na.omit(df$value)))/2)[2],
                       low = '#7bb6d6', na.value = "white") + 
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 10, 
                                   vjust = 1, angle = 45, hjust = 1))




##--------------------------Figure5 D------------------------------------------##

cols = c('#a9479a','#f59771', "#9E6B49", '#48b9b4',
         "#D1DB94", "#87AC3F" ,'#f3c3db', '#86b3e0', "#9E6B49" ,"#E6DEC1" ,"#F3C3DB")

clonal(tregs, 
        reduction = "umap", 
        freq.cutpoint = 2, 
        bins = 10,
        facet = "condition") + 
  scale_color_manual(values = cols)



##--------------------------Supplementary Figure5 D----------------------------##

tregs$cloneType <- factor(tregs$cloneType,
                              levels = c("Single (0 < X <= 1)",
                                         "Small (1 < X <= 5)",
                                         "Medium (5 < X <= 20)", NA))


colorblind_vector <- c("#CEE6F3", "#FFB433", "#FF4B20")


DimPlot(tregs, group.by = "cloneType", order = T) +
  scale_color_manual(values = colorblind_vector, na.value="grey")
