###Analysis for NF2 WT & KO Bulk RNA-seq - Cardiomyocyte DIfferentiation###
require(edgeR)
require(ggplot2)
require(gplots)
require(VennDiagram)
require(RColorBrewer)
require(goseq)
require(scales)
require(reshape)
require(stringr)
require(dplyr)
library(viridis)
require(ggpointdensity)
require(tidyr)
require(ggrepel)

rm(list = ls())

#Set working directory
setwd("C:/Users/lcj_m/OneDrive/Desktop/NF2_KO")

#names(filelist)=c("filename","sample","group")

###Place count files in the working directory - Read Data file###
filelist = data.frame(filename=list.files(pattern="*count.txt"))
filelist$sample = gsub(".count.txt","",filelist$filename)
filelist$group = gsub("_[0-9]","",filelist$sample)
filelist$group = factor(filelist$group,levels=c("WT1_D0","WT2_D0", "KO1_D0", "KO2_D0", "WT1_D1", "WT2_D1", "KO1_D1", "KO2_D1", "WT1_D5", "WT2_D5", "KO1_D5", "KO2_D5", "WT1_D14", "WT2_D14", "KO1_D14", "KO2_D14"))
x <- c("WT1_D0","WT2_D0", "KO1_D0", "KO2_D0", "WT1_D1", "WT2_D1", "KO1_D1", "KO2_D1", "WT1_D5", "WT2_D5", "KO1_D5", "KO2_D5", "WT1_D14", "WT2_D14", "KO1_D14", "KO2_D14")
filelist<-filelist%>%
  arrange(sapply(group, function(y) which(y == x)))

###Run EdgeR - Create count Matrix###
data = readDGE(file=filelist$filename,group=filelist$group,labels=filelist$sample)
filter = apply(data$counts,1,function(x) length(x[x>5])>=2)
filtered = data$counts[ filter, ]
filtered = filtered[grep("^_",rownames(filtered),invert=TRUE),]

###Read Refs file###
gene2biotype = read.table("hg38_gene2biotype.txt")
names(gene2biotype)=c("gid","hgnc","biotype")
genelength = read.table("hg38_genelength.txt")
names(genelength)=c("gid","length")

genelength = genelength[match(rownames(filtered),genelength$gid),]
##This check is necessary to ensure that the gene lengths are accounted for
##Do not delete
table(is.na(genelength$gid))
dim(genelength)
dim(filtered)


grouping = data.frame(sample=filelist$sample,group=filelist$group)
design = model.matrix(~0+group,data=grouping)
colnames(design)=gsub("group","",colnames(design))
design

rnaseq = DGEList(counts=filtered,group=grouping$group)
rnaseq = calcNormFactors(rnaseq)
rnaseq = estimateDisp(rnaseq,design)
fit <- glmFit(rnaseq, design)
normdata = as.data.frame(cpm(rnaseq))
fpkm = as.data.frame(rpkm(rnaseq,gene.length=genelength$length))

###Append HGNC symbol and biotype to ENSG ID###
normdata2 = cbind(normdata,gene2biotype[match(rownames(normdata),gene2biotype$gid),])

###Generate Normalized and FPKM file###
write.table(normdata2,file="normdata_Diff_NF2KO.txt",sep="\t",quote=FALSE)
write.table(fpkm,file="fpkm_diff_NF2KO.txt",sep="\t",quote=FALSE)

contrasts = makeContrasts(D0=(((KO1_D0+KO2_D0)/2 - (WT1_D0+WT2_D0)/2)), D1=(((KO1_D1+KO2_D1)/2 - (WT1_D1+WT2_D1)/2)), D5=(((KO1_D5+KO2_D5)/2 - (WT1_D5+WT2_D5)/2)), D14=(((KO1_D14+KO2_D14)/2 - (WT1_D14+WT2_D14)/2)), levels=design)
contrasts

testDE <- function(x)
{
  test = glmLRT(fit,contrast = contrasts[,x])
  test = as.data.frame(topTags(test,n=100000000))
  test$fpkm.avg = rowMeans(fpkm[match(rownames(test),rownames(fpkm)),])
  test$sig = test$FDR<0.052 & abs(test$logFC)>=1 & test$fpkm>=1
  test = cbind(test,gene2biotype[match(rownames(test),gene2biotype$gid),])
  test = cbind(test,fpkm[match(rownames(test),rownames(fpkm)),])
  return(test)
}


## Generate PLOTs
## Function for Volcano
volcano <- function(x)
{
  x$logP = -log10(x$PValue)
  y = x[ x$sig, ]
  fortext = rbind(head(y[order(y$logFC),],20),tail(y[order(y$logFC),],20))
  
  p <- ggplot(x,aes(x=logFC,y=logP))+geom_point(aes(color=sig))+theme_bw()+theme(panel.grid=element_blank())
  p <- p + scale_color_manual(values=c("grey30","orange"))
  p <- p + xlab("log2FC")+ylab("-log10Pvalue")
  p <- p + geom_text(data=fortext,aes(x=logFC,y=logP,label=hgnc),size=4)
  return(p)
}
## Function for MA
MAplot <- function(x)
{
  p <- ggplot(x,aes(x=logCPM,y=logFC))+geom_point(aes(color=sig))+theme_bw()+theme(panel.grid=element_blank())
  p <- p + scale_color_manual(values=c("grey30","orange"))
  p <- p + xlab("logCPM")+ylab("log2FC")
  return(p)
}

## Function for MA for selected genes ##

GOI<-c("KDR","NKX2-5","ISL1","MEF2C","MESP1", "ETV2", "NANOG","SOX2","PODXL","MYC","CCND1","CTGF","CYR61", "VIM","MSX1","CDX2")
MAplot_genes <- function(x)
{
  fortext <- x %>% filter(hgnc%in%GOI)
  
  p <- ggplot(x,aes(x=logFC,y=logCPM))+geom_point(aes(color=sig))+theme_bw()+theme(panel.grid=element_blank())
  p <- p + xlab("log2FC")+ylab("logCPM")
  p <- p + scale_color_manual(values=c("grey30","red"))+ geom_label_repel(data = fortext, min.segment.length = 0, box.padding = 0.5, max.overlaps= Inf,   
                                                                          aes(label = hgnc))
  return(p)
  
}


## Function for Scatter
scatterplot <- function(x,group1,group2)
{
  temp = data.frame(g1=rowMeans(fpkm[,colnames(fpkm)%in%grouping$sample[grouping$group==group1]]), g2=rowMeans(fpkm[,colnames(fpkm)%in%grouping$sample[grouping$group==group2]]))
  rownames(temp)=rownames(fpkm)
  temp$sig = rownames(temp)%in%rownames(x)[ x$sig ]
  p <- ggplot(temp,aes(x=g1+1,y=g2+1))+geom_point(aes(color=sig))+theme_bw()+theme(panel.grid=element_blank())
  p <- p + scale_color_manual(values=c("grey30","orange"))
  p <- p + xlab(paste("FPKM of",group1))+ylab(paste("FPKM of ",group2))+scale_x_log10()+scale_y_log10()
  return(p)
}

## Function for Heatmap
DEHeat <- function(x,group1,group2,group3,group4)
{
  temp = fpkm[,colnames(fpkm)%in% grouping$sample[grouping$group %in% c(group1,group2,group3,group4)]]
  temp = temp[ rownames(temp) %in% rownames(x)[x$sig],]
  temp = log(temp+1,2)
  temp = (temp-rowMeans(temp))/apply(temp,1,sd)
  plotdata = melt(as.matrix(temp))
  names(plotdata)=c("gid","sample","z")
  plotdata$group = grouping$group[match(plotdata$sample,grouping$sample)]
  plotdata$logFC = x$logFC[match(plotdata$gid,rownames(x))]
  plotdata$genegroup = plotdata$logFC>0
  plotdata$z[ plotdata$z>2]=2
  plotdata$z[ plotdata$z<(-2)]=-2
  p <- ggplot(plotdata,aes(x=sample,y=reorder(gid,logFC)))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank())
  p <- p + scale_fill_gradientn(values=rescale(c(-2,0,2)),colors=c("blue","white","orange"))
  p <- p + theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.spacing=unit(0,"lines")) + facet_grid(genegroup~group,space="free",scale="free")
  return(p)
}




## Function for GOseq
getGeneLists <- function(pwf, goterms, genome, ids){
  gene2cat <- getgo(rownames(pwf), genome, ids)
  cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                    unlist(gene2cat, use.names = FALSE))
  out <- list()
  for(term in goterms){
    tmp <- pwf[cat2gene[[term]],]
    tmp <- rownames(tmp[tmp$DEgenes > 0, ])
    gene2biotype$gid=gsub("\\..*","",gene2biotype$gid)
    out[[term]] <- gene2biotype$hgnc[match(tmp,gene2biotype$gid)]
  }
  out
}


goseqRUN <- function(delist,status,genome)
{
  rownames(delist)=gsub("\\..*","",rownames(delist))
  isSigGene = delist$sig
  if(status == "down"){
    isSigGene[ delist$logFC>0 ] = FALSE
  }else{
    isSigGene[ delist$logFC<0 ] = FALSE
  }
  
  genes = as.integer(isSigGene)
  names(genes)=rownames(delist)
  
  pwf <- nullp(genes,genome,"ensGene")
  GO.wall = goseq(pwf,genome,"ensGene",test.cats=c("GO:BP"))
  ##topGO.wall = head(GO.wall,30)
  topGO.wall = GO.wall[ GO.wall$over_represented_pvalue<0.01, ]
  if(nrow(topGO.wall)>10){
    topGO.wall = head(topGO.wall,10)
  }
  goList = getGeneLists(pwf,topGO.wall$category,genome,"ensGene")
  topGO.wall$EnsemblID <- sapply(topGO.wall$category, function(x) paste0(goList[[x]], collapse = ","))
  return(topGO.wall)
}


plotGOHeatmap <- function(genelist,pathname,plotdata)
{
  tempdata = plotdata[ plotdata$SYMBOL %in% unlist(str_split(genelist,",")), ]
  p1 = ggplot(tempdata,aes(x=sample,y=SYMBOL))+geom_tile(aes(fill=z))+theme_bw()
  p1 = p1 + theme(panel.grid=element_blank(),panel.spacing=unit(0,"lines"),axis.ticks=element_blank())+xlab("")+ylab("")
  p1 = p1 + theme(axis.text.x=element_blank()) + facet_grid(.~group,space="free",scale="free")
  p1 = p1 + scale_fill_gradientn(values=rescale(c(-2,0,2)),colors=c("blue","white","orange"))
  pathname = gsub("\\:","",pathname)
  pathname = gsub(" ","_",pathname)
  filename = paste(pathname,"_GOheat.png",sep="")
  png(filename,width=2500,height=3500,res=300)
  print(p1)
  dev.off()
}

plotGOBAR <- function(GO.wall,direction)
{
  GO.wall$logP = -1*log(GO.wall$over_represented_pvalue,10)
  top10 = head(GO.wall,10)
  p1<-ggplot(top10,aes(x=reorder(term,logP),y=logP))+geom_bar(stat="identity",fill="royalblue")+theme_bw()+theme(panel.grid=element_blank())+coord_flip()+geom_hline(yintercept=seq(0,floor(max(top10$logP)/10)*10,by=10),color="white")+ylab("-log10(P)")+xlab("")
  p1<-p1+ggtitle(direction)
  return(p1)
}

## Function for FGSEA
FGSEA <- function(db,ranks,dbtype){
  ##depending on where you are working with mouse or human
  ##load("pathDB/Mm.c5.symbol.Rdata")
  ##full path to DB
  ## the variable stored is path
  ## we use hgnc remember
  
  pathRdata = paste(dbtype,".RData",sep="")
  pathTxt = paste(dbtype,"_FGSEA.txt",sep="")
  load(db)
  pathgene=c()
  fgseaRes <- fgsea(path,ranks,minSize=15,maxSize=500,nperm=1000)
  sigPATH = fgseaRes[ order(padj, -abs(NES)), ]
  for(i in 1:nrow(sigPATH)){ temp=paste0(unlist(sigPATH[i,8]),collapse=","); pathgene=c(pathgene,temp) }
  
  sigPATH = sigPATH[,1:7]
  sigPATH$gene = pathgene
  return(sigPATH)
}

findSIG <- function(x)
{
  temp = x[ x$sig, ]
  temp$status = as.character(rep("UP"))
  temp$status[ temp$logFC<0 ] = "DOWN"
  return(temp)
}

##This can be used to generate overall heatmap
FinalDEHeat <- function(x)
{
  temp = normdata
  temp = temp[ rownames(temp) %in% x$gid,]
  temp = log(temp+1,2)
  temp = (temp-rowMeans(temp))/apply(temp,1,sd)
  plotdata = melt(as.matrix(temp))
  names(plotdata)=c("gid","sample","z")
  plotdata$group = grouping$group[match(plotdata$sample,grouping$sample)]
  plotdata$genegroup = x$status[match(plotdata$gid,x$gid)]
  
  plotdata$z[ plotdata$z>2]=2
  plotdata$z[ plotdata$z<(-2)]=-2
  p <- ggplot(plotdata,aes(x=sample,y=gid))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank())
  p <- p + scale_fill_gradientn(values=rescale(c(-2,0,2)),colors=c("blue","white","orange"))
  p <- p + theme(axis.text=element_blank(), axis.ticks=element_blank(),panel.spacing=unit(0,"lines")) + facet_grid(genegroup~group,space="free",scale="free")
  return(p)
}


## form plotdata for plotting
require(reshape)
require(ggplot2)
require(scales)
require(ggrepel)

z = (fpkm-rowMeans(fpkm))/apply(fpkm,1,sd)
plotdata = melt(as.matrix(z)) ##formed earlier
names(plotdata) = c("gid","sample","z")
plotdata$group = grouping$group[match(plotdata$sample,grouping$sample)]
plotdata$SYMBOL = gene2biotype$hgnc[match(plotdata$gid,gene2biotype$gid)]
plotdata = plotdata[!is.na(plotdata$SYMBOL),]

####Performing PCA using FPKM values of all groups###
pca = prcomp(t(fpkm),scale=TRUE,center=TRUE)
eigs <- pca$sdev^2
percentage_explained = eigs/sum(eigs)*100
percentage_explained

pcadata = as.data.frame(pca$x)
pcadata$group = filelist$group[match(rownames(pcadata),filelist$sample)]

pcaxl = paste("PC1(",format(round(percentage_explained[1],2),2),"%)",sep="")
pcayl = paste("PC2(",format(round(percentage_explained[2],2),2),"%)",sep="")

pdf("globalpca.pdf")
ggplot(pcadata,aes(x=PC1,y=PC2, label=rownames(pcadata)))+geom_text(vjust=-1.2, size=2)+geom_point(aes(color=group),shape=18,size=3)+theme_bw()+theme(panel.grid=element_blank())+xlab(pcaxl)+ylab(pcayl)+scale_color_manual(values=c("blue","orange", "green", "yellow", "black", "purple", "red", "steelblue", "greenyellow", "deeppink","gray", "lightsalmon", "paleturquoise", "salmon", "turquoise", "violet"))
dev.off()

#Day 1 Comparison KO1/KO2 vs WT1/WT2
D1 = testDE("D1")
write.table(D1,file="D1_DE.txt",sep="\t",quote=FALSE)

D1.heat = DEHeat(D1,"WT1_D1","WT2_D1", "KO1_D1", "KO2_D1")
D1.volcano = volcano(D1)
D1.maplot = MAplot(D1)
D1.scatterplot = scatterplot(D1,"WT2_D1","KO2_D1")
D1.goseq.up = goseqRUN(D1,"up","hg38")
for(i in 1:nrow(D1.goseq.up)){
  print(paste("plotting graph for ",i,"/",nrow(D1.goseq.up)))
  plotGOHeatmap(D1.goseq.up$EnsemblID[i],paste("D1.goseq.up",D1.goseq.up$category[i],D1.goseq.up$term[i]),plotdata)
}
D1.goseq.up.barplot = plotGOBAR(D1.goseq.up,"GOBP UP")
D1.goseq.down = goseqRUN(D1,"down","hg38")
for(i in 1:nrow(D1.goseq.down)){
  print(paste("plotting graph for ",i,"/",nrow(D1.goseq.down)))
  plotGOHeatmap(D1.goseq.down$EnsemblID[i],paste("D1.goseq.down",D1.goseq.down$category[i],D1.goseq.down$term[i]),plotdata)
}
D1.goseq.down.barplot = plotGOBAR(D1.goseq.down,"GOBP DOWN")

pdf("D1.heat_Logfc_1.pdf")
print(D1.heat)
dev.off()

pdf("D1.volcano.pdf")
print(D1.volcano)
dev.off()

pdf("D1.maplot.pdf")
print(D1.maplot)
dev.off()

pdf("D1.scatterplot.pdf")
print(D1.scatterplot)
dev.off()

pdf("D1.goseq.up.pdf")
print(D1.goseq.up.barplot)
dev.off()

pdf("D1.goseq.down.pdf")
print(D1.goseq.down.barplot)
dev.off()


#Plot overall heatmap with hiererical clustering 

filelist$color[ filelist$group=="WT1_D0" ] = "red"
filelist$color[ filelist$group=="WT2_D0" ] = "brown"
filelist$color[ filelist$group=="KO1_D0" ] = "blue"
filelist$color[ filelist$group=="KO2_D0"] = "turquoise"
filelist$color[ filelist$group=="WT1_D1"] = "red"
filelist$color[ filelist$group=="WT2_D1"] = "brown"
filelist$color[ filelist$group=="KO1_D1"] = "blue"
filelist$color[ filelist$group=="KO2_D1"] = "turquoise"
filelist$color[ filelist$group=="WT1_D5"] = "red"
filelist$color[ filelist$group=="WT2_D5"] = "brown"
filelist$color[ filelist$group=="KO1_D5"] = "blue"
filelist$color[ filelist$group=="KO2_D5"] = "turquoise"
filelist$color[ filelist$group=="WT1_D14"] = "red"
filelist$color[ filelist$group=="WT2_D14"] = "brown"
filelist$color[ filelist$group=="KO1_D14"] = "blue"
filelist$color[ filelist$group=="KO2_D14"] = "turquoise"



##To plot heatmap for specific genes - Extended Figure 6B ##

fpkm_wgenelist= cbind(fpkm,gene2biotype[match(rownames(fpkm),gene2biotype$gid),])
fpkm2<-fpkm_wgenelist
PL<- c("NANOG", "SOX2", "POU5F1", "ZFP42", "PODXL")
ME<- c("GSC","EOMES","T","MIXL1","FOXA2","SOX17","BMP4","CXCR4","PDGFRA","EVX1","LHX1","FOXC1","TBR1","SNAI1", "GATA4", "GATA6", "HNF4A", "FOXA1", "SOX7")
CPC<- c("MESP1","HAND1","TBX1","ETV2","ISL1","KDR","PDGFRA","TWIST1")
NE<- c("GBX2","OTX2","SOX1","PAX6","NEUROG3","NEFM","SGCA","BMP4","SCGE","FABP3")
CM<- c("TBX5", "NPPA", "MYH7", "ACTN2", "MYL7", "NKX2-5", "ACTA1", "MYH6", "MYL2", "ACTC1", "MEF2C", "TNNT2","GATA4")
YAP<-c("NPPB","MYC", "CTGF", "CYR61", "VIM", "CCND1")
GOI<-c("KDR","NKX2-5","ISL1","MEF2C","MESP1", "ETV2", "NANOG","SOX2","PODXL","MYC","CCND1","CTGF","CYR61", "VIM","MSX1","CDX2")

filtered_PL<-fpkm2 %>% filter(hgnc %in% PL)
filtered_ME<-fpkm2 %>% filter(hgnc %in% ME)
filtered_CPC<-fpkm2 %>% filter(hgnc %in% CPC)
filtered_CM<-fpkm2 %>% filter(hgnc %in% CM)
filtered_NE<-fpkm2 %>% filter(hgnc %in% NE)
filtered_YAP<-fpkm2 %>% filter(hgnc %in% YAP)

PL_Filtered<-select(filtered_PL, -c(gid, hgnc,biotype))
ME_Filtered<-select(filtered_ME, -c(gid, hgnc, biotype))
CPC_Filtered<-select(filtered_CPC, -c(gid, hgnc, biotype))
CM_Filtered<-select(filtered_CM, -c(gid, hgnc, biotype))
NE_Filtered<-select(filtered_NE, -c(gid, hgnc, biotype))
YAP_Filtered<-select(filtered_YAP, -c(gid, hgnc, biotype))

write.table(fpkm_wgenelist2,file="fpkm_wgenelist_DE.txt",sep="\t",quote=FALSE)
write.table(fpkm2,file="fpkm2_frame.txt",sep="\t",quote=FALSE, row.names=FALSE)

#### To plot NE-specific genes #####
fpkm_NE<- NE_Filtered
heatmapdegene = paste("degene_heatmap_NE_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=1500,height=2300,res=300)
fpkm_NE = (fpkm_NE-rowMeans(fpkm_NE))/apply(fpkm_NE,1,sd)
rownames(fpkm_NE) = D1$hgnc[match(rownames(fpkm_NE),rownames(D1))]
rownames(fpkm_CM) = make.names(DoxGroup$hgnc[match(rownames(DoxGroup_fpkm.sig),rownames(DoxGroup))], unique=TRUE)
hc = hclust(as.dist(1-cor(fpkm_NE,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(fpkm_NE),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(fpkm_NE),Rowv=as.dendrogram(hr), Colv=FALSE,scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","orange"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color), margins=c(8,8))
dev.off()

#### To plot PL-specific genes #####
fpkm_PL<- PL_Filtered
heatmapdegene = paste("degene_heatmap_PL_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=4000,height=1500,res=300)
fpkm_PL = log(fpkm_PL+1,2)
fpkm_PL = (fpkm_PL-rowMeans(fpkm_PL))/apply(fpkm_PL,1,sd)
rownames(fpkm_PL) = D1$hgnc[match(rownames(fpkm_PL),rownames(D1))]
#rownames(fpkm_ME) = make.names(DoxGroup$hgnc[match(rownames(DoxGroup_fpkm.sig),rownames(DoxGroup))], unique=TRUE)
hc = hclust(as.dist(1-cor(fpkm_PL,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(fpkm_PL),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(fpkm_PL),Rowv=as.dendrogram(hr), Colv=FALSE, scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","red"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color), margins=c(8,8), rowsep=1:nrow(fpkm_PL), colsep=c(6,12, 17, 23, 28, 34, 39), sepcolor="black")
dev.off()


#### To plot CM-specific genes #####
fpkm_CM<- CM_Filtered
fpkm_CM= log(fpkm_CM+1,2)
heatmapdegene = paste("degene_heatmap_CM_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=3000,height=1500,res=300)
fpkm_CM = (fpkm_CM-rowMeans(fpkm_CM))/apply(fpkm_CM,1,sd)
rownames(fpkm_CM) = D1$hgnc[match(rownames(fpkm_CM),rownames(D1))]
#rownames(fpkm_ME) = make.names(DoxGroup$hgnc[match(rownames(DoxGroup_fpkm.sig),rownames(DoxGroup))], unique=TRUE)
hc = hclust(as.dist(1-cor(fpkm_CM,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(fpkm_CM),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(fpkm_CM),Rowv=as.dendrogram(hr), Colv=FALSE, scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","red"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color), margins=c(8,8), rowsep=1:nrow(fpkm_CM), colsep=c(6,12, 17, 23, 28, 34, 39), sepcolor="black")
dev.off()

#### To plot YAP-target genes #####
fpkm_YAP<- YAP_Filtered
fpkm_YAP = log(fpkm_YAP+1,2)
heatmapdegene = paste("degene_heatmap_YAP_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=4000,height=1500,res=300)
fpkm_YAP = (fpkm_YAP-rowMeans(fpkm_YAP))/apply(fpkm_YAP,1,sd)
rownames(fpkm_YAP) = D1$hgnc[match(rownames(fpkm_YAP),rownames(D1))]
#rownames(fpkm_ME) = make.names(DoxGroup$hgnc[match(rownames(DoxGroup_fpkm.sig),rownames(DoxGroup))], unique=TRUE)
hc = hclust(as.dist(1-cor(fpkm_YAP,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(fpkm_YAP),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(fpkm_YAP),Rowv=as.dendrogram(hr), Colv=FALSE, scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","red"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color), margins=c(8,8), rowsep=1:nrow(fpkm_YAP), colsep=c(6,12, 17, 23, 28, 34, 39), sepcolor="black")
dev.off()

#### To plot ME-specific genes #####
fpkm_ME<- ME_Filtered
fpkm_ME= log(fpkm_ME+1,2)
heatmapdegene = paste("degene_heatmap_ME_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=3000,height=2000,res=300)
fpkm_ME = (fpkm_ME-rowMeans(fpkm_ME))/apply(fpkm_ME,1,sd)
rownames(fpkm_ME) = D1$hgnc[match(rownames(fpkm_ME),rownames(D1))]
#rownames(fpkm_ME) = make.names(DoxGroup$hgnc[match(rownames(DoxGroup_fpkm.sig),rownames(DoxGroup))], unique=TRUE)
hc = hclust(as.dist(1-cor(fpkm_ME,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(fpkm_ME),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(fpkm_ME),Rowv=as.dendrogram(hr), Colv=FALSE, scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","red"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color), margins=c(8,8), rowsep=1:nrow(fpkm_ME), colsep=c(6,12, 17, 23, 28, 34, 39), sepcolor="black")
dev.off()


#### To plot CPC-specific genes #####
fpkm_CPC<-CPC_Filtered
fpkm_CPC= log(fpkm_CPC+1,2)
heatmapdegene = paste("degene_heatmap_CPC_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=3000,height=1500,res=300)
fpkm_CPC = (fpkm_CPC-rowMeans(fpkm_CPC))/apply(fpkm_CPC,1,sd)
rownames(fpkm_CPC) = D1$hgnc[match(rownames(fpkm_CPC),rownames(D1))]
#rownames(fpkm_ME) = make.names(DoxGroup$hgnc[match(rownames(DoxGroup_fpkm.sig),rownames(DoxGroup))], unique=TRUE)
hc = hclust(as.dist(1-cor(fpkm_CPC,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(fpkm_CPC),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(fpkm_CPC),Rowv=as.dendrogram(hr), Colv=FALSE, scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","red"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color), margins=c(8,8), rowsep=1:nrow(fpkm_CM), colsep=c(6,12, 17, 23, 28, 34, 39), sepcolor="black")
dev.off()

#### To plot all DE genes with reference to Day 1 #####
fpkm_D1<- fpkm
heatmapdegene = paste("degene_heatmap_D1_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=1500,height=2300,res=300)
D1_fpkm.sig = fpkm_D1[ rownames(fpkm_D1) %in% rownames(D1)[ D1$sig ], ]
D1_fpkm.sig = (D1_fpkm.sig-rowMeans(D1_fpkm.sig))/apply(D1_fpkm.sig,1,sd)
rownames(D1_fpkm.sig) = D1$hgnc[match(rownames(D1_fpkm.sig),rownames(D1))]
hc = hclust(as.dist(1-cor(D1_fpkm.sig,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(D1_fpkm.sig),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(D1_fpkm.sig), Rowv=as.dendrogram(hr), Colv=FALSE,scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","orange"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color))
dev.off()


#Day 5 Comparison KO1/KO2 vs WT1/WT2
D5 = testDE("D5")
write.table(D5,file="D5_DE_logFC_1.txt",sep="\t",quote=FALSE)

D5.heat = DEHeat(D5,"WT1_D5","WT2_D5", "KO1_D5", "KO2_D5")
D5.volcano = volcano(D5)
D5.maplot = MAplot(D5)
D5.scatterplot = scatterplot(D5,"WT2_D5","KO2_D5")
D5.goseq.up = goseqRUN(D5,"up","hg38")
for(i in 1:nrow(D5.goseq.up)){
  print(paste("plotting graph for ",i,"/",nrow(D5.goseq.up)))
  plotGOHeatmap(D5.goseq.up$EnsemblID[i],paste("D5.goseq.up",D5.goseq.up$category[i],D5.goseq.up$term[i]),plotdata)
}
D5.goseq.up.barplot = plotGOBAR(D5.goseq.up,"GOBP UP")
D5.goseq.down = goseqRUN(D5,"down","hg38")
for(i in 1:nrow(D5.goseq.down)){
  print(paste("plotting graph for ",i,"/",nrow(D5.goseq.down)))
  plotGOHeatmap(D5.goseq.down$EnsemblID[i],paste("D5.goseq.down",D5.goseq.down$category[i],D5.goseq.down$term[i]),plotdata)
}
D5.goseq.down.barplot = plotGOBAR(D5.goseq.down,"GOBP DOWN")

pdf("D5.heat.pdf")
print(D5.heat)
dev.off()

pdf("D5.volcano.pdf")
print(D5.volcano)
dev.off()

pdf("D5.maplot.pdf")
print(D5.maplot)
dev.off()

pdf("D5.scatterplot.pdf")
print(D5.scatterplot)
dev.off()

pdf("D5.goseq.up.pdf")
print(D5.goseq.up.barplot)
dev.off()

pdf("D5.goseq.down.pdf")
print(D5.goseq.down.barplot)
dev.off()


###Plot heatmap with hiererical clustering for Day 5 ###

fpkm_D5<- fpkm
heatmapdegene = paste("degene_heatmap_D5_combined_FDR_0.05",".png",sep="")
png(heatmapdegene,width=1500,height=2300,res=300)
D5_fpkm.sig = fpkm_D5[ rownames(fpkm_D5) %in% rownames(D5)[ D5$sig ], ]
D5_fpkm.sig = (D5_fpkm.sig-rowMeans(D5_fpkm.sig))/apply(D5_fpkm.sig,1,sd)
rownames(D5_fpkm.sig) = make.names(D5$hgnc[match(rownames(D5_fpkm.sig),rownames(D5))], unique=TRUE) 
hc = hclust(as.dist(1-cor(D5_fpkm.sig,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(D5_fpkm.sig),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(D5_fpkm.sig),Colv=as.dendrogram(hc),Rowv=as.dendrogram(hr),scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","orange"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color))
dev.off()


##Day 14 Comparison KO1/KO2 vs WT1/WT2
D14 = testDE("D14")
write.table(D14,file="D14_DE.txt",sep="\t",quote=FALSE)

D14.heat = DEHeat(D14,"WT1_D14","WT2_D14", "KO1_D14", "KO2_D14")
D14.volcano = volcano(D14)
D14.maplot = MAplot(D14)
D14.scatterplot = scatterplot(D14,"WT2_D14","KO2_D14")
D14.goseq.up = goseqRUN(D14,"up","hg38")
for(i in 1:nrow(D14.goseq.up)){
  print(paste("plotting graph for ",i,"/",nrow(D14.goseq.up)))
  plotGOHeatmap(D14.goseq.up$EnsemblID[i],paste("D14.goseq.up",D14.goseq.up$category[i],D14.goseq.up$term[i]),plotdata)
}
D14.goseq.up.barplot = plotGOBAR(D14.goseq.up,"GOBP UP")
D14.goseq.down = goseqRUN(D14,"down","hg38")
for(i in 1:nrow(D14.goseq.down)){
  print(paste("plotting graph for ",i,"/",nrow(D14.goseq.down)))
  plotGOHeatmap(D14.goseq.down$EnsemblID[i],paste("D14.goseq.down",D14.goseq.down$category[i],D14.goseq.down$term[i]),plotdata)
}
D14.goseq.down.barplot = plotGOBAR(D14.goseq.down,"GOBP DOWN")


pdf("D14.heat_.pdf")
print(D14.heat)
dev.off()

pdf("D14_.volcano.pdf")
print(D14.volcano)
dev.off()

pdf("D14.maplot.pdf")
print(D14.maplot)
dev.off()

pdf("D14.scatterplot.pdf")
print(D14.scatterplot)
dev.off()

pdf("D14.goseq.up.pdf")
print(D14.goseq.up.barplot)
dev.off()

pdf("D14.goseq.down.pdf")
print(D14.goseq.down.barplot)
dev.off()


#Plot heatmap with hiererical clustering 

fpkm_D14<- fpkm
heatmapdegene = paste("degene_heatmap_D14_FDR_0.05",".png",sep="")
png(heatmapdegene,width=1500,height=2300,res=300)
D14_fpkm.sig = fpkm_D14[ rownames(fpkm_D14) %in% rownames(D14)[ D14$sig ], ]
D14_fpkm.sig = (D14_fpkm.sig-rowMeans(D14_fpkm.sig))/apply(D14_fpkm.sig,1,sd)
rownames(D14_fpkm.sig) = make.names(D14$hgnc[match(rownames(D14_fpkm.sig),rownames(D14))], unique=TRUE)
hc = hclust(as.dist(1-cor(D14_fpkm.sig,method="spearman")),method="ward.D2")
hr = hclust(as.dist(1-cor(t(D14_fpkm.sig),method="pearson")),method="ward.D2")
heatmap.2(as.matrix(D14_fpkm.sig),Colv=as.dendrogram(hc),Rowv=as.dendrogram(hr),scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","orange"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color))
dev.off()


##Find significant DE genes from D1, D5 and D14##
D1.sig=findSIG(D1)
D5.sig=findSIG(D5)
D14.sig=findSIG(D14)



###Produce merged Heatmap
union.list = unique(c(as.character(rownames(D1.sig)),as.character(rownames(D5.sig)), as.character(rownames(D14.sig))))
final.frame = data.frame(gid=union.list)

#Set1 comparing D1
#Set2 comparing D5
#Set3 comparing D14
final.frame$Set1 = D1.sig$status[match(final.frame$gid,rownames(D1.sig))]
final.frame$Set2 = D5.sig$status[match(final.frame$gid,rownames(D5.sig))]
final.frame$Set3 = D14.sig$status[match(final.frame$gid,rownames(D14.sig))]
final.frame$status = paste(final.frame$Set1,final.frame$Set2, final.frame$Set3)
#setting up final.frame
final.frame$hgnc = D1$hgnc[match(final.frame$gid, rownames(D1))]
final.frame$logFC.D1 = D1$logFC[match(final.frame$gid, rownames(D1))]
final.frame$sig.D1 = D1$sig[match(final.frame$gid, rownames(D1))]

final.frame$logFC.D5 = D5$logFC[match(final.frame$gid, rownames(D5))]
final.frame$sig.D5 = D5$sig[match(final.frame$gid, rownames(D5))]

final.frame$logFC.D14 = D14$logFC[match(final.frame$gid, rownames(D14))]
final.frame$sig.D14 = D14$sig[match(final.frame$gid, rownames(D14))]
table(final.frame$status)

write.table(table(final.frame$status), file="table_changes.txt", sep="\t", quote=FALSE)

DOWN_DOWN_DOWN<-filter(final.frame, status=="DOWN DOWN DOWN")
DOWN_DOWN_NA<-filter(final.frame, status=="DOWN DOWN NA")
DOWN_NA_DOWN<-filter(final.frame, status=="DOWN NA DOWN")
DOWN_NA_NA<-filter(final.frame, status=="DOWN NA NA") 
DOWN_UP_UP<-filter(final.frame, status=="DOWN UP UP")
NA_DOWN_DOWN<-filter(final.frame, status=="NA DOWN DOWN")
NA_NA_DOWN<-filter(final.frame, status=="NA NA DOWN") 
NA_DOWN_NA<-filter(final.frame, status=="NA DOWN NA")
NA_NA_UP<-filter(final.frame, status=="NA NA UP")
NA_UP_NA<-filter(final.frame, status=="NA UP NA")
NA_UP_UP<-filter(final.frame, status=="NA UP UP")
UP_NA_NA<-filter(final.frame, status=="UP NA NA")
UP_NA_UP<-filter(final.frame, status=="UP NA UP")
UP_UP_NA<-filter(final.frame, status=="UP UP NA")
UP_UP_UP<-filter(final.frame, status=="UP UP UP")



write.table(final.frame,file="DE_Changes.txt",sep="\t",quote=FALSE)
write.table(DOWN_DOWN_DOWN,file="DE_DOWN_DOWN_DOWN.txt",sep="\t",quote=FALSE)
write.table(DOWN_DOWN_NA,file="DE_DOWN_DOWN_NA.txt",sep="\t",quote=FALSE)
write.table(DOWN_NA_DOWN,file="DE_DOWN_NA_DOWN.txt",sep="\t",quote=FALSE)
write.table(DOWN_NA_NA,file="DE_DOWN_NA_NA.txt",sep="\t",quote=FALSE)
write.table(DOWN_UP_UP,file="DE_DOWN_UP_UP.txt",sep="\t",quote=FALSE)
write.table(NA_DOWN_DOWN,file="DE_NA_DOWN_DOWN.txt",sep="\t",quote=FALSE)
write.table(NA_NA_DOWN,file="DE_NA_NA_DOWN.txt",sep="\t",quote=FALSE)
write.table(NA_DOWN_NA,file="DE_NA_DOWN_NA.txt",sep="\t",quote=FALSE)
write.table(NA_NA_UP,file="DE_NA_NA_UP.txt",sep="\t",quote=FALSE)
write.table(NA_UP_NA,file="DE_NA_UP_NA.txt",sep="\t",quote=FALSE)
write.table(NA_UP_UP,file="DE_NA_UP_UP.txt",sep="\t",quote=FALSE)
write.table(UP_NA_NA,file="DE_UP_NA_NA.txt",sep="\t",quote=FALSE)
write.table(UP_NA_UP,file="DE_UP_NA_UP.txt",sep="\t",quote=FALSE)
write.table(UP_UP_NA,file="DE_UP_UP_NA.txt",sep="\t",quote=FALSE)
write.table(UP_UP_UP, file="DE_UP_UP_UP.txt", sep="\t", quote=FALSE)


##
pdf("merged_comparison.pdf",width=5,height=5)
print(FinalDEHeat(final.frame))
dev.off()

###plot special maplot
ggmaplot(D1, Xval = D1$logCPM,
         Yval = D1$logFC,
         basic.Color = "black",
         title.Text = "Day 1 - Mesendoderm Induction",
         Ylim = c(-6,6),
         use.PointDensity = TRUE)


####Prepare dataframe for Triage analysis####

fpkm2<-select(fpkm2, -c(biotype,gid))
fpkm2<-fpkm2 %>% select(hgnc, everything())
rownames(fpkm2)<-NULL
names(fpkm2)[names(fpkm2) == "hgnc"] <- "Gene"

# Remove duplicates based on Gene
fpkm2<-fpkm2[!duplicated(fpkm2$Gene), ]
fpkm_D1_WT<-select(fpkm2,c("Gene","WT1_D1_1","WT1_D1_2", "WT2_D1_1", "WT2_D1_2", "WT2_D1_3"))
fpkm_D1_KO<-select(fpkm2,c("Gene","KO1_D1_1","KO1_D1_2", "KO1_D1_3", "KO2_D1_1", "KO2_D1_2", "KO2_D1_3"))

fpkm_D5_WT<-select(fpkm2,c("Gene","WT1_D5_1","WT1_D5_2", "WT2_D5_1", "WT2_D5_2", "WT2_D5_3"))
fpkm_D5_KO<-select(fpkm2,c("Gene","KO1_D5_1","KO1_D5_2", "KO1_D5_3", "KO2_D5_1", "KO2_D5_2", "KO2_D5_3"))

fpkm_D14_WT<-select(fpkm2,c("Gene","WT1_D14_1","WT1_D14_2", "WT2_D14_1", "WT2_D14_2", "WT2_D14_3"))
fpkm_D14_KO<-select(fpkm2,c("Gene","KO1_D14_1","KO1_D14_2", "KO1_D14_3", "KO2_D14_1", "KO2_D14_2", "KO2_D14_3"))


###CSV files for TRIAGE###
write.csv(fpkm_D1_WT, file="fpkm_D1_WT.csv", row.names = FALSE)
write.csv(fpkm_D1_KO, file="fpkm_D1_KO.csv", row.names = FALSE)
write.csv(fpkm_D5_KO, file="fpkm_D5_KO.csv", row.names = FALSE)
write.csv(fpkm_D5_WT, file="fpkm_D5_WT.csv", row.names = FALSE)
write.csv(fpkm_D14_WT, file="fpkm_D14_WT.csv", row.names = FALSE)
write.csv(fpkm_D14_KO, file="fpkm_D14_KO.csv", row.names = FALSE)



###END###