############## FigS7C
#obtain CG site in certain genomic regions using bedtools

input <- file.path("./")
files <- list.files(input,"txt$")
library(dplyr)
table <- c()
for (i in files){
    print(i)
    sample <- gsub(".cov5-.+txt","",i)
    sample <- gsub("Foxp3cre-cfp1-KOA_1_val_1_bismark_bt2_pe.deduplicated.bismark","KO",sample)
    sample <- gsub("Foxp3cre-cfp1-WTA_1_val_1_bismark_bt2_pe.deduplicated.bismark","WT",sample)
    repeats <- gsub(".+cov5-mm9_","",i)
    repeats <- gsub(".txt","",repeats)
    #print(repeats)
    #print(sample)
    x <-read.table(i)
    colnames(x) <- c("chr","start","end","meth","coverage")
    x <- x %>% mutate(treat=sample,repeats=repeats)
    table <- rbind(table,x)
}

library(ggplot2)
pdf("meth.pdf",width=10,height=4) 
table$treat <- factor(table$treat,levels=c("WT","KO"))
ggplot(table,aes(x=repeats,y=meth,fill=treat)) + geom_violin(scale="width") + theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1)) + xlab("") + scale_fill_manual(values=c("blue","red")) +
 stat_summary(fun=mean, geom="point", position=position_dodge(0.9), pch=4,color="black", size=3) + ylab("mCG/CG") 
dev.off()

