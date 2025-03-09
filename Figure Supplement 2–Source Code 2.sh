# FigS2C


computeMatrix reference-point -S Foxp3cre-WT-IgG_treat_pileup.bw Foxp3cre-WT-CFP1_treat_pileup.bw -R CGI.bed -b 5000 -a 5000 -o CGI.CFP1.gz --numberOfProcessors 30 --referencePoint center

plotProfile -m CGI.CFP1.gz -out CGI.CFP1.profile.pdf --samplesLabel IgG CFP1 --regionsLabel CGI --perGroup

