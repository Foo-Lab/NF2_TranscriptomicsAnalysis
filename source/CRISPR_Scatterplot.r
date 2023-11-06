## Script-Plot Scatterplot for CRISPR Screen - Figure 1B ##

pacman::p_load(tidyverse, readxl, janitor, 
               ggfortify, ggrepel,gridExtra,
               pheatmap)
pacman::p_unload(plyr)

rm(list = ls())

setwd("C:/Users/lcj_m/OneDrive/Desktop/CRISPRScreen")

zscore <- read.table("zscore2.txt")
zscore_pvalue <-read_excel("zscore_pvalue.txt")
filtered<-zscore %>% filter(`GFP Negative_1`>1 & `GFP Negative_2`>1 & `GFP Positive_1`>1 & `GFP Positive_2`>1)


significant  <- zscore_pvalue %>% filter(pvalue>1 & zscore>1)
significant <- as.data.frame(significant)
nonsignificant <- zscore_pvalue %>% filter(!pvalue>1 & zscore>1)
  
zscore_pvalue<-as.data.frame(zscore_pvalue)
control<-zscore_pvalue %>% filter(gene=="Control")
control<-as.data.frame(control)
genetgfb<-c("ACVR1","ACVR1B","BAMBI")
tgfb<- zscore_pvalue %>% filter(gene%in%genetgfb)
genetgfyap<-c("BMPR1A", "BMPR2", "SMAD2")
tgfyap<- zscore_pvalue %>% filter(gene%in%genetgfyap)
yap<-zscore_pvalue %>% filter(gene=="FBXW11")
Wnt<-zscore_pvalue %>% filter(gene=="GSK3B")
NF2<-zscore_pvalue %>% filter(gene=="NF2")

label<-c("USP34", "CARM1", "AP2M1", "AP2S1", "GSK3B", "TBX20", "ACVR1", "ACVR1B", "BMPR1A", "BMPR2", "SMAD2", "MIXL1", "ARID1A", "NF2", "JMJD1C")
tolabel<-zscore_pvalue %>% filter(gene%in%label)
pdf("Rplot_Label_final.pdf", width=8, height=6)
ggplot(zscore_pvalue, aes(x=zscore, y=pvalue))+
  geom_point(col="black", alpha=0.25, size=1)+
  geom_abline(intercept=0, slope=0, col="black")+
  geom_abline(intercept=1, slope=0, col="blue", lty=2, lwd=1)+
  geom_point(data=control, aes(col="Control"), alpha=0.25, size=1)+
  geom_point(data=nonsignificant,aes(col="Non-enriched Genes"), alpha= 0.25, size =1)+
  geom_point(data=significant, aes(col="Positive Enrichment GFP-Negative"), alpha=0.8, size=2)+
  scale_color_manual(name=NULL, values=c("Positive Enrichment GFP-Negative"="red", "Control"="orange", "Borderline Hits"="Blue", "Non-enriched Genes"="grey"))+
  theme_bw()+
  labs(x = "Z Score", y = " -LOG10(P)")+
  geom_label_repel(data = tolabel, min.segment.length = 0, box.padding = unit(0.5, "lines"), label.size=0.1, fill = alpha(c("white"),0.5),
                   aes(label = gene))+
  geom_label_repel(data = borderline2, min.segment.length = 0, box.padding = unit(0.5, "lines"),  label.size=0.1, fill = alpha(c("white"),0.5),
                   aes(label = gene))+
  guides(colour = guide_legend(override.aes = list(size=2)))


dev.off()
