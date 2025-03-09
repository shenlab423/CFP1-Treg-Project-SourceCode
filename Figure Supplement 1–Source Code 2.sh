
##### FigS1B

computeMatrix reference-point -S Foxp3cre-WT-H3K4me3_treat_pileup.bw -R H3K4me3_peaks.bed -b 5000 -a 5000 -o H3K4me3.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m H3K4me3.gz -out H3K4me3.heatmap.pdf --samplesLabel Foxp3-WT-H3K4me3 --regionsLabel H3K4me3_peaks --colorList "white,blue" --yMin 0 --missingDataColor white


##### FigS1C


computeMatrix reference-point -S Tconv_H3K4me3.bw Treg_H3K4me3.bw -R Treg_specific_promoter.bed -b 5000 -a 5000 -o Tregspcific_H3K4me3.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m Tregspcific_H3K4me3.gz -out Tregspcific_H3K4me3.heatmap.pdf --samplesLabel Tconv Treg --regionsLabel Treg_specific_promoter --colorList "white,blue" --yMin 0 --missingDataColor white


