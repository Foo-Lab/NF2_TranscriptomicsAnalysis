setwd("C:/Users/lcj_m/OneDrive/Desktop/NF2_KO")
library(ggplot2)
library(dplyr)
library(ggrepel)


###Figure 5B-C TRIAGE visualization###
###Input WT Discordance and expression data, and label TF and HK genes###
DS = read.table("DS.txt", header=TRUE, row.names=NULL)
TF = read.table("TF.txt", header=TRUE,row.names=NULL )
HK = read.table("HK_genes.txt", header=TRUE,row.names=NULL )
RTS = read.table("rts.txt", header=TRUE,row.names=NULL )
EXP = read.table("wtexpression_fpkm.txt", header=TRUE,row.names=NULL)

colnames(DS)=c("gene","DS")
#Rank DS SCORE
DS$TRIAGE_RANK[order(-DS$DS)] <-1:nrow(DS)
write.table(DS,file="DS_RANK_WT.txt",sep="\t",quote=FALSE)

DS$tf = TF$type[match(DS$gene,TF$gene)]
DS$hk = HK$type[match(DS$gene,HK$gene)]
DS$RTS = RTS$RTS[match(DS$gene, RTS$gene)]
DS$EXP_RANK = EXP$Rank[match(DS$gene, EXP$Gene)]

###SUPPLEMENTARY TABLE 5####
EXP1 = read.table("EXPRANK.txt", header=TRUE,row.names=NULL)
TRIAGE1 = read.table("TRIAGERANK.txt", header=TRUE, row.names=NULL)
Supp5<-left_join(TRIAGE1,EXP1, 'Gene')
write.table(Supp5,file="Supplementary Table 5.txt",sep="\t",quote=FALSE)

DS$type = paste(DS$tf,DS$hk)
table(DS$type)


DS<-distinct(DS,gene, .keep_all= TRUE)
TF<-DS %>% filter(type=="TF NA")
HK<-DS %>% filter(type=="NA HK")
tolabel<- DS%>%filter(DS>0.5)
tolabel<-tolabel%>%filter(tolabel$TRIAGE_RANK<13)

p<-ggplot(DS, aes(x=EXP_RANK, y=DS))+
  geom_point(col="black", alpha=0.25, size=1)+
  geom_point(data=TF, aes(col="Regulatory Genes"), alpha=1, size=2)+
  geom_point(data=HK,aes(col="Housekeeping Genes"), alpha= 1, size =2)+
  scale_color_manual(name=NULL, values=c("Regulatory Genes"= "red", "Housekeeping Genes"="blue"))+
  geom_label_repel(data = tolabel, min.segment.length = 0, box.padding = 0.5, max.overlaps= Inf,aes(label = gene))+
  theme_bw()+theme(panel.grid=element_blank())+
  labs(x = "Gene Expression Ranking", y = "Discordance Score")
ggsave(filename="DS_plot_expRANK_WTD14.pdf", plot=p, width=10, height=5, units="in")


q<-ggplot(HK, aes(x=gene, y=DS))+
  geom_point(alpha=1, size=1.5,aes(col="Housekeeping Genes"))+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_color_manual(name=NULL, values=c("Housekeeping Genes"= "blue"))+
  labs(x = "genes", y = "Discordance Score")+ylim(0,3)
ggsave(filename="DS_plot_HK.pdf", plot=q, width=10, height=5, units="in")

