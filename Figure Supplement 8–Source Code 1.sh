#### FigS8A


computeMatrix reference-point -S Foxp3cre-WT-Foxp3_treat_pileup.bw Foxp3cre-KO-Foxp3_treat_pileup.bw -R GSE121279_peaks-foxp3.mm9.bed -b 5000 -a 5000 -o Foxp3binding.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m Foxp3binding.gz -out Foxp3binding.heatmap.pdf --samplesLabel Foxp3-WT Foxp3-KO --regionsLabel Foxp3peaks --colorList "white,blue" --yMin 0 --missingDataColor white




