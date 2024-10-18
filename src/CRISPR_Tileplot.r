pacman::p_load(tidyverse, readxl, janitor, 
               ggfortify, ggrepel,gridExtra,
               pheatmap)
pacman::p_unload(plyr)
require(ggside)
require(ggplot2)
require(RColorBrewer)

rm(list = ls())

###Set your working directory here where you save the folder##
setwd("C:/Users/Mick Lee/Desktop/CRISPRTileplot")

zscore <- read.table("zscore2.txt",header=TRUE)
zscore <- as.data.frame(zscore)
z <- zscore


##########################################################################

###Plotting the distribution of your library, setting X scale consistent to -3 to 3###

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

overall

dev.off()

#########################################################################

###############PLOT non-targetting Controls in orange#######################

GOI_control="Control"
NTC<-z %>% filter(gene%in%GOI_control)
pdf("zscore_density_plot_Controls.pdf")
overall<-ggplot(NTC, aes(zscore)) +
  geom_tile( # geom_tile() draws rectangular colored areas
    aes(
      y = 1, # draw all tiles at the same y location
      fill = after_stat(density)  # use computed density for fill
    ),
    stat = "density",
    n = 200,
    show.legend = F# number of points calculated by stat_density() 
  )+scale_fill_gradient(low="white",high="orange")+ scale_y_continuous(expand = c(0,0))+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),panel.background = element_blank(), plot.background = element_blank())+xlim(c(-3,3))+xlab("")+ylab("")+theme(ggside.panel.border = element_blank())

overall
dev.off()


### Tileplot of selected candidate gene of interest - Figure 1C ###
##Example for ARID1A and JMJD1C##


GOI<-c("ARID1A", "JMJD1C")
list<-z %>% filter(gene%in%GOI)

######Example for to plot Selected genes ###############
##Make sure you keep the x-scale consistent


GOI<-c("ARID1A", "JMJD1C")
list<-z %>% filter(gene%in%GOI)

##Export plot in png ##
plot = paste("zscore_enriched_ARID1A_JMJD1C","_",".png",sep="")
png(plot,width=2000,height=500,res=300)
p<-test+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background = element_blank(),
                         panel.border = element_rect(colour = "black", fill=NA, size=1),
                         plot.background = element_blank())+xlim(c(-3,3))+geom_vline(data = list, aes(xintercept=as.numeric(zscore)), color ="red", lwd=1, lty=1)+xlab("")+ylab("")+coord_equal()+facet_grid(list$gene~., switch="y")
p





dev.off()

