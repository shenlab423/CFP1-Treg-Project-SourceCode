### Fig1D
#H3K4me3

computeMatrix reference-point -S Foxp3cre-WT-H3K4me3.bw -R HM_clustered_promoter.nTreg.bed -b 5000 -a 5000 -o H3K4me3.cluster.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m H3K4me3.cluster.gz -out H3K4me3.cluster.heatmap.pdf --samplesLabel H3K4me3 --regionsLabel Foxp3peaks --colorList "white,blue" --yMin 0 --missingDataColor white

#H3K27me3

computeMatrix reference-point -S Foxp3cre-WT-H3K27me3.bw -R HM_clustered_promoter.nTreg.bed -b 5000 -a 5000 -o H3K27me3.cluster.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m H3K27me3.cluster.gz -out H3K27me3.cluster.heatmap.pdf --samplesLabel H3K27me3 --regionsLabel Foxp3peaks --colorList "white,blue" --yMin 0 --missingDataColor white


#FOXP3

computeMatrix reference-point -S Treg-FOXP3.bw -R HM_clustered_promoter.nTreg.bed -b 5000 -a 5000 -o Foxp3.cluster.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m Foxp3.cluster.gz -out Foxp3.cluster.heatmap.pdf --samplesLabel Foxp3 --regionsLabel Foxp3peaks --colorList "white,blue" --yMin 0 --missingDataColor white


