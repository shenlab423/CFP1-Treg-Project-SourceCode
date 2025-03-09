########## FigS7G


computeMatrix reference-point -S Foxp3cre-WT-H3K4me1.bw Foxp3cre-KO-H3K4me1.bw -R H3K4me1_peaks.bed -b 5000 -a 5000 -o H3K4me1.peaks.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m H3K4me1.peaks.gz -out H3K4me1.peaks.heatmap.pdf --samplesLabel Foxp3-WT Foxp3-KO --regionsLabel H3K4me1_peaks --colorList "white,blue" --yMin 0 --missingDataColor white



########## FigS7H-I
###### obtain top 5% longest H3K4me3 peaks as broad domains
###### obtain broad domains covered genes using bedtools

computeMatrix reference-point -S Foxp3cre-WT-H3K4me3.bw Foxp3cre-KO-H3K4me3.bw -R H3K4me3_broad_domain_genes.bed -b 5000 -a 5000 -o H3K4me3.BD.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m H3K4me3.BD.gz -out H3K4me3.BD.heatmap.pdf --samplesLabel Foxp3-WT Foxp3-KO --regionsLabel H3K4me3_BD --colorList "white,blue" --yMin 0 --missingDataColor white



computeMatrix reference-point -S Foxp3cre-IgG.bw Foxp3cre-WT-CFP1.bw -R H3K4me3_broad_domain_genes.bed -b 5000 -a 5000 -o CFP1.BD.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m CFP1.BD.gz -out CFP1.BD.heatmap.pdf --samplesLabel IgG CFP1 --regionsLabel H3K4me3_BD --colorList "white,blue" --yMin 0 --missingDataColor white

