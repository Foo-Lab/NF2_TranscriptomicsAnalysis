## Tradeseq - Pseudotime Analysis ##
## No. of replicates per group required to be consistent ##
require(edgeR)
require(gplots)
require(VennDiagram)
require(RColorBrewer)
require(goseq)
require(scales)
require(reshape)
require(stringr)
require(dplyr)
require(circlize)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(clusterExperiment)
library(ComplexHeatmap)
library(GetOptLong)

rm(list =ls())

setwd("C:/Users/lcj_m/OneDrive/Desktop/NF2_KO/tradeseq")

#names(filelist)=c("filename","sample","group")

##Setting the Matrix input for tradeseq##
filelist = data.frame(filename=list.files(pattern="*count.txt"))
filelist$sample = gsub(".count.txt","",filelist$filename)
filelist$group = gsub("_[0-9]","",filelist$sample)

x <- c("WT1_D0", "WT2_D0", "WT1_D1", "WT2_D1", "WT1_D5", "WT2_D5", "WT1_D14", "WT2_D14", "KO1_D0", "KO2_D0", "KO1_D1", "KO2_D1", "KO1_D5", "KO2_D5", "KO1_D14", "KO2_D14")
filelist<-filelist%>%
  arrange(sapply(group, function(y) which(y == x)))

data = readDGE(file=filelist$filename,group=filelist$group,labels=filelist$sample)
filter = apply(data$counts,1,function(x) length(x[x>5])>=2)
filtered = data$counts[ filter, ]
filtered = filtered[grep("^_",rownames(filtered),invert=TRUE),]
temp<-as.data.frame(filtered)
rownames(temp)= gsub("\\..*","",rownames(temp))
colnames(temp)=gsub("_D1_","_1_",colnames(temp))
colnames(temp)=gsub("_D0_","_0_",colnames(temp))
colnames(temp)=gsub("_D5_","_5_",colnames(temp))
colnames(temp)=gsub("_D14_","_14_",colnames(temp))
matrix<-as.matrix(temp)


grouping = data.frame(sample=filelist$sample,group=filelist$group)
design = model.matrix(~0+group,data=grouping)
colnames(design)=gsub("group","",colnames(design))
design


## Data input from Refs ##
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
##


rnaseq = DGEList(counts=filtered,group=grouping$group)
rnaseq = calcNormFactors(rnaseq)
rnaseq = estimateDisp(rnaseq,design)
fit <- glmFit(rnaseq, design)
normdata = as.data.frame(cpm(rnaseq))
fpkm = as.data.frame(rpkm(rnaseq,gene.length=genelength$length))
rownames(fpkm)= gsub("\\..*","",rownames(fpkm))

## Building the count matrix ##
counts<- matrix
time <- as.numeric(unlist(lapply(strsplit(colnames(counts),split="_"),"[",2)))
treat <- factor(unlist(lapply(strsplit(colnames(counts),split="_"),"[",1)), levels=c("WT", "KO"))

# filter counts: at least 5 counts in 3 replicats
keep <- rowSums(counts>5)>=3
counts <- counts[keep,]
counts<-na.omit(counts)

## Tradeseq - Trajectory inference analysis ## WT VS KO ##
library(tradeSeq)
time <- matrix(time, nrow=ncol(counts), ncol=2, byrow=FALSE)
rownames(time) <- colnames(counts)
weights <- matrix(0, nrow=ncol(counts), ncol=2)
weights[1:(ncol(counts)/2), 1] <- 1
weights[(ncol(counts)/2+1):ncol(counts), 2] <- 1

## evaluate optimal K
infMat <- evaluateK(counts, pseudotime=time, cellWeights=weights, nGenes=250, k=3:4, verbose = T)

## fit GAM
gamList <- fitGAM(counts, pseudotime=time, cellWeights=weights, nknots=4)

## look at fit for six random genes
set.seed(81)
id <- sample(1:length(gamList), size=6)
par(mfrow=c(2,3))
for(ii in 1:6) plotSmoothers(gamList[[id[ii]]], main=resPat[id[ii],"pvalue"], ylim=c(3,9))


time2<-as.data.frame(time)

#### RUN Slingshot
set.seed(7)
pseudotime <- slingPseudotime(time2, na = FALSE)
cellWeights <- slingCurveWeights(time)

## Pattern test for different gene expression pattern ##
resPat <- patternTest(gamList, nPoints=8)
resPat2_wpvalue <- patternTest(gamList, nPoints=8)

## Adding statistical significance to pattern test ##
hist(resPat$pvalue)
sum(p.adjust(resPat$pvalue, "fdr") <= 0.01, na.rm=TRUE)

# Comparing Berge et al., Nature Comms. 2020, For neuronal differentiation of DE genes, Pairwise comparisons were made with DESeq2 and found 7486 DE genes. We have 7184, which is very comparable.
deTradeSeq <- rownames(counts)[p.adjust(resPat$pvalue, "fdr") <= 0.01]
resPat2_wpvalue$FDR<- p.adjust(resPat2_wpvalue$pvalue, "fdr") <= 0.05
patternde$FDR<-p.adjust(resPat2_wpvalue$pvalue, "fdr") <= 0.05
patternde$FDR_2<-p.adjust(resPat2_wpvalue$pvalue, "fdr") <= 0.01

## plot six most significant genes.
oo <- order(resPat$pvalue, decreasing=FALSE)
par(mfrow=c(2,3))
for(ii in seq(1,6)) plotSmoothers(gamList[[oo[ii]]], main=rownames(resPat)[oo[ii]], ylim=c(1,9))
dev.new()
par(mfrow=c(2,3))
for(ii in seq(7,12)) plotSmoothers(gamList[[oo[ii]]], main=rownames(resPat)[oo[ii]], ylim=c(1,9))

## cluster

deGenes <- rownames(counts)[p.adjust(resPat$pvalue, "fdr") <= 0.01]
deGenes <- deGenes[deGenes %in% names(gamList)]
clusRes2 <- clusterExpressionPatterns(gamList, nPoints=8, genes=deGenes, nReducedDims=4)
degenes<-as.data.frame(deGenes)
degenes$gid = degenes$deGenes
degenes = select(degenes, -1)

##Match gene2biotype to degenes derived from pattern test, Supplementary Table 6##
gene2biotype = read.table("hg38_gene2biotype.txt")
names(gene2biotype)=c("gid","hgnc","biotype")
gene2biotype$gid=gsub("\\..*","",gene2biotype$gid)
degenes = cbind(degenes,gene2biotype[match(degenes$gid,gene2biotype$gid),])
write.table(degenes,file="tradeseq_degenes_patterntest.txt",sep="\t",quote=FALSE)

patternde = cbind(resPat2_wpvalue,gene2biotype[match(rownames(resPat2_wpvalue),gene2biotype$gid),])
write.table(patternde,file="tradeseq_patterntest.txt",sep="\t",quote=FALSE)

clusterLabels <- primaryCluster(clusRes2$rsec)

cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
7
##Plotting of clusters - Figure 5D and Extended Figure 7D ##

for (xx in cUniq[11]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:8,
                                y = rep(range(clusRes2$yhatScaled[cId, ]),
                                        8 / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusRes2$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:8, 2),
                                  y = clusRes2$yhatScaled[geneId, ],
                                  lineage = rep(0:1, each = 8)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("orange", "darkseagreen3"),
                       breaks = c("0", "1"))  
  print(p)
}


cId <- which(clusterLabels == 8)
p <- ggplot(data = data.frame(x = 1:8,
                              y = rep(range(clusRes2$yhatScaled[cId, ]),
                                      8 / 2)),
            aes(x = x, y = y)) +
  geom_point(alpha = 0) +
  labs(title = paste0("Cluster ", 35),  x = "Pseudotime", y = "Normalized expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))


geneId <- rownames(clusRes2$yhatScaled)[cId[35]]
p <- p +
  geom_line(data = data.frame(x = rep(1:8, 2),
                              y = clusRes2$yhatScaled[geneId, ],
                              lineage = rep(0:1, each = 8)),
            aes(col = as.character(lineage), group = lineage), lwd = 1.5)

p <- p + guides(color = FALSE) +
  scale_color_manual(values = c("orange", "darkseagreen3"),
                     breaks = c("0", "1"))  


print(p)

## Plotting of Smoothers for CDX2 ENSG00000165556 ## 

yhat <- predictCells(models = gamList, gene = "ENSG00000165556")

ysmooth <- predictSmooth(models = gamList, gene = "ENSG00000165556", nPoints = 8)



## Plotting of Smoothers MSX1 ENSG00000163132 ##

ysmooth <- predictSmooth(models = gamList, gene = "ENSG00000163132", nPoints = 8)

## Plotting of Smoothers NKX2-5 ENSG00000183072 ##

ysmooth <- predictSmooth(models = gamList, gene = "ENSG00000183072", nPoints = 8)

## Plotting of Smoothers SOX2 ENSG00000181449 ##

ysmooth <- predictSmooth(models = gamList, gene = "ENSG00000181449", nPoints = 8)

## Association test - Supplementary Table 7 ##

assocRes <- associationTest(gamList, lineages = TRUE)
assocRes = cbind(assocRes,gene2biotype[match(rownames(assocRes),gene2biotype$gid),])

write.table(assocRes,file="tradeseq_assoctest.txt",sep="\t",quote=FALSE)

## Pattern test with wald statistics ##
library(ggplot2)

resPat$Gene <- rownames(resPat)
resPat$pattern <- resPat$waldStat
resPat <- resPat[, c("Gene", "pattern")]

resPat$Gene <- rownames(resPat)
resPat$end <- resPat$waldStat
resPat <- resPat[, c("Gene", "end")]

compare <- merge(patternRes, resPat, by = "Gene", all = FALSE)
compare$transientScore <- 
  rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2

ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "patternTest Wald Statistic (log scale)",
       y = "diffEndTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()

defaultMar<-par("mar")
plotCMar<-c(1.1,8.1,4.1,1.1)
par(mar=plotCMar)
plotClusters(clusRes2$rsec,main="Clusters from RSEC", whichClusters="workflow", colData=c("Biological_Condition","Published2"), axisLine=-1)


##plot heat map of top 5000 associated genes - Figure 5E##
##Data input from Pseudotime folder##
### Read only WT FPKM which will be used as anchor ###

fpkm_WT = read.table("fpkm_WT.txt") #fpkm only with wildtype group

top5k = read.table("top5k_Assoc.txt")
colnames(top5k) = c("gid", "hgnc")
top5k<-as.data.frame(top5k)

#Plot heatmap with hiererical clustering 
filelist$color = as.character(rep("blue"))
filelist$color[ filelist$group=="WT2_D0" ] = "red"
filelist$color[ filelist$group=="WT2_D1" ] = "orange"
filelist$color[ filelist$group=="WT2_D5"] = "yellow"
filelist$color[ filelist$group=="WT1_D0" ] = "red"
filelist$color[ filelist$group=="WT1_D1" ] = "orange"
filelist$color[ filelist$group=="WT1_D5"] = "yellow"
filelist$color[ filelist$group=="KO1_D0" ] = "red"
filelist$color[ filelist$group=="KO1_D1" ] = "orange"
filelist$color[ filelist$group=="KO1_D5"] = "yellow"
filelist$color[ filelist$group=="KO2_D0" ] = "red"
filelist$color[ filelist$group=="KO2_D1" ] = "orange"
filelist$color[ filelist$group=="KO2_D5"] = "yellow"
heatmapdegene = paste("AssocHeatmap_WT_KO_5k",".png",sep="")
png(heatmapdegene,width=1500,height=2300,res=300)
fpkm = (fpkm-rowMeans(fpkm))/apply(fpkm,1,sd)
rownames(fpkm) = make.names(top5k$hgnc[match(rownames(fpkm),top5k$gid)], unique=TRUE)
rownames(fpkm_WT) = make.names(top5k$hgnc[match(rownames(fpkm_WT),top5k$gid)], unique=TRUE) #aim to cluster by WT first
fpkm <- fpkm[grepl("^NA", rownames(fpkm))==F,]
fpkm_WT <- fpkm_WT[grepl("^NA", rownames(fpkm_WT))==F,]

#fpkm_WT <- fpkm[ -c(20:41) ]
hc = hclust(as.dist(1-cor(fpkm,method="spearman")),method="ward.D")
hr = hclust(as.dist(1-cor(t(fpkm_WT),method="pearson")),method="ward.D")
pdf("AssocHeatmap_top5k.pdf")
heatmap.2(as.matrix(fpkm),Rowv=as.dendrogram(hr), Colv=FALSE,scale="none",breaks=c(seq(-2,2,by=0.04)),col=colorRampPalette(c("blue","white","orange"))(100),density.info="none",trace="none",ColSideColors=as.character(filelist$color), margins=c(8,8))
dev.off()

###try plot in complex heatmap###
# Your list of genes
genelist <- c("TNNT2", "TBX1", "T", "SOX2", "SNAI1", "POU5F1", "MYH7", "MYH6", "MSX1", "ISL1", "FOXA2", "NODAL", "FOXC1", "NPPA", "NANOG", "CXCR4", "PDGFRA", "EOMES", "TWIST1", "HAND2", "HAND1", "MYH6", "MYH7" ) 
# Labels to plot
as.matrix(fpkm)
labels <- rownames(fpkm) 
# Set the rest of genes to ""
labels[!labels %in% genelist] <- "" 
my_colors = c("blue", "white", "orange")
col_fun<-colorRamp2(c(-2,0,2), c("blue","white", "orange"))
my_colors = colorRampPalette(my_colors)(25) 

colnames(filelist) <- c("filename", "sample", "group", "color")
ann <- data.frame(filelist$group)
colnames(ann)<-c("group")
colours <- list("group" = c("WT1_D0" = "red","WT1_D14" = "blue", "WT1_D5"="yellow", "WT1_D1"="orange",
                            "WT2_D0" = "red","WT2_D14" = "blue", "WT2_D5"="yellow", "WT2_D1"="orange",
                            "KO1_D0" = "red","KO1_D14" = "blue", "KO1_D5"="yellow", "KO1_D1"="orange",
                            "KO2_D0" = "red","KO2_D14" = "blue", "KO2_D5"="yellow", "KO2_D1"="orange"))

colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

column_split = rep("Wild-Type", 41)
column_split[20:41] = "NF2 Knockout"


pdf("heatmap_5kassoc_orange.pdf")
pheatmap(as.matrix(fpkm),
         scale="row",
         color = col_fun,
         cluster_rows = hr,
         cluster_cols = FALSE, 
         column_split = column_split,
         top_annotation=colAnn,
         cutree_rows = 4) +
  rowAnnotation(mark = anno_mark(labels = test, at=c(1445,	1571,	192,	1155,	3044,	86,	4226,	555	,2741,	3628,	1919,	4368,	1518,	3825,	3076,	3942,	1330,	1462), which="row"))
dev.off()
test<-c("CXCR4"	,"FOXA2","FOXC1","HAND1","HAND2",	"ISL1",	"MYH6",	"MYH7",	"NODAL"	,"NPPA",	"PDGFRA",	"POU5F1",	"SNAI1",	"SOX2",	"T",	"TBX1",	"TNNT2",	"TWIST1")

write.table(test,file="labels.txt",sep="\t",quote=FALSE)


##END##