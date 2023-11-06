pacman::p_load(tidyverse, readxl, janitor, 
               ggfortify, ggrepel,gridExtra,
               pheatmap)
pacman::p_unload(plyr)
require(ggside)
require(ggplot2)
require(RColorBrewer)

rm(list = ls())

setwd("C:/Users/lcj_m/OneDrive/Desktop/CRISPRScreen")

zscore <- read.table("zscore2.txt")
zscore <- as.data.frame(zscore)
z <- zscore


### Tileplot of selected candidate gene of interest - Figure 1C ###
##Example for ARID1A and JMJD1C##

z$gene2=z$gene
GOI<-c("ARID1A", "JMJD1C")
list<-z %>% filter(gene%in%GOI)
#list$gene2 = list$gene


test<-ggplot(z, aes(zscore)) +
  geom_tile( # geom_tile() draws rectangular colored areas
    aes(
      y = 1, # draw all tiles at the same y location
      fill = after_stat(density)  # use computed density for fill
    ),
    stat = "density",
    n = 200,
    show.legend = F# number of points calculated by stat_density() 
  )+scale_fill_gradient(low="white",high="black")+ scale_y_continuous(expand = c(0,0))
p<-test+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),
                         panel.border = element_rect(colour = "black", fill=NA, size=1),
                         plot.background = element_blank())+xlim(c(-3,3))+geom_vline(data = list, aes(xintercept=as.numeric(zscore)), color ="red", lwd=1, lty=1)+xlab("")+ylab("")+coord_equal()+facet_grid(list$gene2~., switch="y")

plot = paste("zscore_enriched_ARID1A_JMJD1C","_",".png",sep="")
png(plot,width=2000,height=500,res=300)
pdf("zscore_ARID1A_JMJD1C.pdf")
p

dev.off()

####PLOT non-targetting Controls###
overall<-ggplot(list, aes(zscore)) +
  geom_tile( # geom_tile() draws rectangular colored areas
    aes(
      y = 1, # draw all tiles at the same y location
      fill = after_stat(density)  # use computed density for fill
    ),
    stat = "density",
    n = 200,
    show.legend = F# number of points calculated by stat_density() 
  )+scale_fill_gradient(low="white",high="orange")+ scale_y_continuous(expand = c(0,0))+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),panel.background = element_blank(), plot.background = element_blank())+xlim(c(-3,3))+xlab("")+ylab("")+theme(ggside.panel.border = element_blank())
pdf("zscore_density_plot_Controls.pdf")
overall
dev.off()


plot = paste("zscore_density_plot","_",".png",sep="")
png(plot,width=2000,height=500,res=300)
overall<-ggplot(z, aes(zscore)) +
  geom_tile( # geom_tile() draws rectangular colored areas
    aes(
      y = 1, # draw all tiles at the same y location
      fill = after_stat(density)  # use computed density for fill
    ),
    stat = "density",
    n = 200,
    show.legend = F# number of points calculated by stat_density() 
  )+scale_fill_gradient(low="white",high="black")+ scale_y_continuous(expand = c(0,0))+geom_xsidedensity(aes(y = after_stat(density)),show.legend = F)+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background = element_blank(), plot.background = element_blank())+xlim(c(-3,3))+xlab("")+ylab("")+theme(ggside.panel.border = element_blank())
pdf("zscore_density_plot.pdf")
overall
dev.off()

