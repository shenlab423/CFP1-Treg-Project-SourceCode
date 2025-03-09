###### Fig7A-B
###### TSSs were grouped into broad TSS(covered by >5kb H3K4me3 peaks), medium TSS(covered by 1-5 kb peaks),narrow TSS(covered by < 1kb peaks)

computeMatrix reference-point -S Foxp3cre-WT-H3K4me3.bw Foxp3cre-KO-H3K4me3.bw -R H3K4me3_broad_TSS.bed H3K4me3_medium_TSS.bed H3K4me3_narrow_TSS.bed -b 5000 -a 5000 -o H3K4me3.TSSgroups.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m H3K4me3.TSSgroups.gz -out H3K4me3.TSSgroups.heatmap.pdf --samplesLabel Foxp3-WT Foxp3-KO --regionsLabel broadTSS mediumTSS narrowTSS --colorList "white,blue" --yMin 0 --missingDataColor white




computeMatrix reference-point -S Foxp3cre-IgG.bw Foxp3cre-WT-CFP1.bw -R H3K4me3_broad_TSS.bed H3K4me3_medium_TSS.bed H3K4me3_narrow_TSS.bed -b 5000 -a 5000 -o CFP1.TSSgroups.gz --numberOfProcessors 30 --referencePoint center

plotHeatmap -m CFP1.TSSgroups.gz -out CFP1.TSSgroups.heatmap.pdf --samplesLabel IgG CFP1 --regionsLabel broadTSS mediumTSS narrowTSS --colorList "white,blue" --yMin 0 --missingDataColor white

