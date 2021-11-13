
setwd("")
library(GEOquery)
library(dplyr)
library(limma)
library(ComplexHeatmap)
library(RColorBrewer)
library(clusterProfiler)
library(tibble)
library(ggplot2)
library(cowplot)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
input=read.table("input.txt",header = T,sep = "\t")

input <- as.matrix(input)
input[is.na(input)] <- "Unknown"
input <- as.data.frame(input)
input$riskScore <- as.numeric(input$riskScore)
input$Histology <- factor(input$Histology,levels = c("astrocytoma","glioblastoma","oligoastrocytoma","oligodendroglioma"))
input$IDH.status <- factor(input$IDH.status,levels = c("Mutant","Unknown","WT"))
input$Gender <- factor(input$Gender,levels = c("female","male"))
input$X1p.19q.codeletion <- factor(input$X1p.19q.codeletion,levels = c("codel","non-codel","Unknown"))
input$IDH.codel.subtype <- factor(input$IDH.codel.subtype,levels = c("IDHmut-codel","IDHmut-non-codel","IDHwt","Unknown"))
input$Grade <- factor(input$Grade,levels = c("G2","G3","G4"))
input$index <- 1:nrow(input)


darkred   <- "#F2042C"
blue      <- "#0A0745"
lightgrey <- "#dcddde"

My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0,0.1,0,0.1, "cm"))

p1 <- ggplot(input, aes(index, riskScore))+
  geom_bar(stat="identity",col = blue)+
  My_Theme1 +
  labs(y = "riskScore") +
  scale_x_continuous(expand = c(0,0)) #不留空

p2 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = Histology))+
  My_Theme2+
  labs(y = "Histology")+
  scale_fill_manual(values = c("#74C065","#B096C2","#FAA96C","#EFF882")) +
  scale_x_continuous(expand = c(0,0)) 

p3 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = IDH.status))+
  My_Theme2+
  labs(y = "IDH.status")+
  scale_fill_manual(values = c("#84CBA8","#F2FA9C","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) 

p4 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = Gender))+
  My_Theme2+
  labs(y = "Gender")+
  scale_fill_manual(values=c("#E00F0A","#3D6197")) +
  scale_x_continuous(expand = c(0,0)) 

p5 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = X1p.19q.codeletion))+
  My_Theme2 +
  labs(y = "X1p.19q.codeletion")+
  scale_fill_manual(values=c("#64B685","#FC6D4C","#6D6466")) +
  scale_x_continuous(expand = c(0,0)) 

p6 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = IDH.codel.subtype))+
  My_Theme2 +
  labs(y = "IDH.codel.subtype")+
  scale_fill_manual(values=c("#C0FF3E","#BFEFFF","#BC8F8F","#6D6466")) +
  scale_x_continuous(expand = c(0,0))

p7 <- ggplot(input,aes(index,1))+
  geom_tile(aes(fill = Grade))+
  My_Theme2 +
  labs(y = "Grade")+
  scale_fill_manual(values=c("#D2691E","#CDCD00","#EE2C2C")) +
  scale_x_continuous(expand = c(0,0)) 

legend_a <- get_legend(p2+theme(legend.position = "bottom"))
legend_b <- get_legend(p3+theme(legend.position = "bottom"))
legend_c <- get_legend(p4+theme(legend.position = "bottom"))
legend_d <- get_legend(p5+theme(legend.position = "bottom"))
legend_e <- get_legend(p6+theme(legend.position = "bottom"))
legend_f <- get_legend(p7+theme(legend.position = "bottom"))
p <- plot_grid(p1,p2,p3,p4,p5,p6,p7,
               align = "v",axis = "l",
               ncol = 1, rel_heights = c(4,1,1,1,1),
               legend_a,legend_b,legend_c,legend_d,legend_e,legend_f)
p

