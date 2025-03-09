######## Fig2A

computeMatrix scale-regions -S IgG_treat_pileup.bw CFP1-chip_treat_pileup.bw -R mm9_TSS.bed -m 5000 -b 5000 -a 5000 -o CFP1_TSS.gz --numberOfProcessors 30 --missingDataAsZero

plotProfile -m CFP1_TSS.gz -out CFP1_TSS.profile.pdf --samplesLabel IgG CFP1 --regionsLabel mm9.TSS --perGroup


######## Fig2C

computeMatrix reference-point -S Foxp3-chip_treat_pileup.bw CFP1-chip_treat_pileup.bw -R Foxp3-specific.bed Foxp3-CFP1-overlap.bed CFP1-specific.bed -b 5000 -a 5000 -o Peaks-Foxp3-CFP1.gz --numberOfProcessors 30 --referencePoint center --missingDataAsZero

plotHeatmap -m Peaks-Foxp3-CFP1.gz -out Peaks-Foxp3-CFP1.pdf --samplesLabel Foxp3 CFP1 --regionsLabel Foxp3-specific.peaks Foxp3-CFP1-overlap.peaks CFP1-specific.peaks --colorList "white,red"


