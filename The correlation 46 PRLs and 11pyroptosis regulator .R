

library(reshape2)
library(RColorBrewer)

options(stringsAsFactors = F)

setwd("")        
up <- read.table("pvalue.txt",sep = "\t",check.names = F,header = T,row.names=1)  
dn <- read.table("cor.txt",sep = "\t",check.names = F,header = T,row.names=1)     
dn=t(dn)
up=t(up)

colVector=c("#AB221F","#3878C1","#FFFADD")


gene.level <- as.character(rownames(dn)) 
cancer.level <- as.character(colnames(dn))

dn.long <- setNames(melt(dn), c('Gene', 'Cancer', 'Frequency'))
dn.long$Categrory <- "DN"
up.long <- setNames(melt(up), c('Gene', 'Cancer', 'Frequency'))
up.long$Categrory <- "UP"

dn.long$range <- cut(dn.long$Frequency, 
                     breaks = seq(floor(min(dn.long$Frequency)),
                                  ceiling(max(dn.long$Frequency)),0.01))
rangeMat1 <- levels(dn.long$range) 
rbPal1 <- colorRampPalette(colors = c(colVector[3],"white",colVector[1]))
col.vec1 <- rbPal1(length(rangeMat1)); names(col.vec1) <- rangeMat1
dn.long$color <- col.vec1[as.character(dn.long$range)]


up.long$range <- cut(up.long$Frequency, breaks = seq(floor(min(up.long$Frequency)),ceiling(max(up.long$Frequency)),0.01)) 
rangeMat2 <- levels(up.long$range)
rbPal2 <- colorRampPalette(colors = c(colVector[3],colVector[2]))
col.vec2 <- rbPal2(length(rangeMat2)); names(col.vec2) <- rangeMat2
up.long$color <- col.vec2[as.character(up.long$range)]


heatmat <- rbind.data.frame(dn.long,up.long) 
pdf("heatmap2.pdf",width = 9,height = 4)
layout(mat=matrix(c(1,0,1,2,1,0,1,3,1,0),5,2,byrow=T),widths=c(length(cancer.level),2))


par(bty="n", mgp = c(2,0.5,0), mar = c(5.1, 5.5, 3, 3),tcl=-.25,xpd = T)
x=as.numeric(factor(heatmat$Cancer,levels = cancer.level))
y=as.numeric(factor(heatmat$Gene,levels = gene.level))

plot(1,xlim=c(1,length(unique(x))+1),ylim=c(1,length(unique(y))+1),
     xaxs="i", yaxs="i",xaxt="n",yaxt="n",
     type="n",bty="n",xlab="",ylab="",
     main = "Coexpression across cancer types",cex.main=2)

for(i in 1:nrow(heatmat)) {
    if(heatmat$Categrory[i]=="DN") polygon(x[i]+c(0,1,1),y[i]+c(0,0,1),col=heatmat$color[i]) #填充左上角
    if(heatmat$Categrory[i]=="UP") {
        polygon(x[i]+c(0,1,0),y[i]+c(0,1,1),col=heatmat$color[i]) #填充右下角
        if(heatmat$Frequency[i]<0.001){
            text(x[i]+0.5,y[i]+0.8,"***",cex=0.8)
        }else if(heatmat$Frequency[i]<0.01){
            text(x[i]+0.5,y[i]+0.8,"**",cex=0.8)
        }else if(heatmat$Frequency[i]<0.05){
            text(x[i]+0.5,y[i]+0.8,"*",cex=0.8)
        }
    }
}

axis(1,at = sort(unique(x)) + 0.5,labels = cancer.level,lty = 0,las = 2)  #添加x轴坐标
axis(2,at = sort(unique(y)) + 0.5,labels = gene.level,lty = 0,las = 1)    #添加y轴坐标
mtext("Cancer types",side = 1,line = 3.5,cex=1.2)    #x轴名称


par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
barplot(rep(1,length(col.vec2)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec2)),horiz=TRUE,
        axes = F, col=col.vec2)  # Loss
axis(4,at=c(1,ceiling(length(col.vec2)/2),length(col.vec2)),c(round(min(up),1),'Pvalue',round(max(up),1)),tick=FALSE)
par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
barplot(rep(1,length(col.vec1)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec1)),horiz=TRUE,
	axes = F, col=col.vec1)  # Gain
axis(4,at=c(1,ceiling(length(col.vec1)/2),length(col.vec1)),c(round(min(dn),1),'Cor',round(max(dn),1)),tick=FALSE)
dev.off()


