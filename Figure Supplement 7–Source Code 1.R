################# FigS7A

#calculate average methylation levels around TSS-TES using Source Code 2

x1 <- read.table("Foxp3cre-cfp1-WTA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.TSS.upstream.txt")
x2 <- read.table("Foxp3cre-cfp1-WTA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.TSS.txt")
x3 <- read.table("Foxp3cre-cfp1-WTA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.TSS.downstream.txt")


library(dplyr)
x2 <- x2 %>% mutate(start=ifelse(V6=="+",V2,(V3+10000)),end=ifelse(V6=="+",V3,(V2-10000)))

x <- merge(x1,x2,by.x = c("V1","V3","V4","V5","V6"),by.y = c("V1","start","V4","V5","V6"))

x <- merge(x,x3, by.x = c("V1","end","V4","V5","V6"),by.y = c("V1","V2","V4","V5","V6"))

x <- x[,c(8:57,61:160,163:212)]

meanrmna <- function(x){
  return(mean(x,na.rm = T))
}

WT <- apply(x, 2, meanrmna)



x1 <- read.table("Foxp3cre-cfp1-KOA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.TSS.upstream.txt")
x2 <- read.table("Foxp3cre-cfp1-KOA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.TSS.txt")
x3 <- read.table("Foxp3cre-cfp1-KOA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.TSS.downstream.txt")


library(dplyr)
x2 <- x2 %>% mutate(start=ifelse(V6=="+",V2,(V3+10000)),end=ifelse(V6=="+",V3,(V2-10000)))

x <- merge(x1,x2,by.x = c("V1","V3","V4","V5","V6"),by.y = c("V1","start","V4","V5","V6"))

x <- merge(x,x3, by.x = c("V1","end","V4","V5","V6"),by.y = c("V1","V2","V4","V5","V6"))

x <- x[,c(8:57,61:160,163:212)]

meanrmna <- function(x){
  return(mean(x,na.rm = T))
}

KO <- apply(x, 2, meanrmna)


WT <- data.frame(x=1:200,y=WT,type="WT")
KO <- data.frame(x=1:200,y=KO,type="KO")

methylation <- rbind(WT,KO)

library(ggplot2)
methylation$type <- factor(methylation$type,levels = c("WT","KO"))
ggplot(methylation,aes(x=x,y=y,col=type)) + geom_line() + theme_classic() + xlab("") + ylab("mCG/CG (%)") + scale_color_manual(values = c("blue","red"))


################# FigS7B
#calculate average methylation levels of 10kb bins using Source Code 2

x1 <- read.table("Foxp3cre-cfp1-WTA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.10kb.bed")
x2 <- read.table("Foxp3cre-cfp1-KOA_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov5.10kb.bed")

x3 <- merge(x1,x2,by=c("V1","V2","V3"))

library(LSD)
heatscatter(x3$V4.x,x3$V4.y,xlab="WT",ylab = "CFP1 KO",main="10kb bins")

