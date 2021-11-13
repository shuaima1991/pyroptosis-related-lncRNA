

library(reshape2)
library(RColorBrewer)

options(stringsAsFactors = F)

setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/46lncRNA与11焦亡基因相关性")        #设置工作目录
up <- read.table("pvalue.txt",sep = "\t",check.names = F,header = T,row.names=1)  #读取左上角的数据
dn <- read.table("cor.txt",sep = "\t",check.names = F,header = T,row.names=1)     #读取右下角数据
dn=t(dn)
up=t(up)
#设置颜色
colVector=c("#AB221F","#3878C1","#FFFADD")

#行名和列名
gene.level <- as.character(rownames(dn)) 
cancer.level <- as.character(colnames(dn))

#把行转为列
dn.long <- setNames(melt(dn), c('Gene', 'Cancer', 'Frequency'))
dn.long$Categrory <- "DN"
up.long <- setNames(melt(up), c('Gene', 'Cancer', 'Frequency'))
up.long$Categrory <- "UP"

#右下角颜色
dn.long$range <- cut(dn.long$Frequency, 
                     breaks = seq(floor(min(dn.long$Frequency)),
                                  ceiling(max(dn.long$Frequency)),0.01))
rangeMat1 <- levels(dn.long$range) # 提出分割区间
rbPal1 <- colorRampPalette(colors = c(colVector[3],"white",colVector[1]))
col.vec1 <- rbPal1(length(rangeMat1)); names(col.vec1) <- rangeMat1
dn.long$color <- col.vec1[as.character(dn.long$range)]

#左上角颜色
up.long$range <- cut(up.long$Frequency, breaks = seq(floor(min(up.long$Frequency)),ceiling(max(up.long$Frequency)),0.01)) 
rangeMat2 <- levels(up.long$range)
rbPal2 <- colorRampPalette(colors = c(colVector[3],colVector[2]))
col.vec2 <- rbPal2(length(rangeMat2)); names(col.vec2) <- rangeMat2
up.long$color <- col.vec2[as.character(up.long$range)]

#合并右下角和左上角
heatmat <- rbind.data.frame(dn.long,up.long) #汇总热图矩阵
pdf("heatmap2.pdf",width = 9,height = 4)
layout(mat=matrix(c(1,0,1,2,1,0,1,3,1,0),5,2,byrow=T),widths=c(length(cancer.level),2))

#热图绘制区域
par(bty="n", mgp = c(2,0.5,0), mar = c(5.1, 5.5, 3, 3),tcl=-.25,xpd = T)
x=as.numeric(factor(heatmat$Cancer,levels = cancer.level))
y=as.numeric(factor(heatmat$Gene,levels = gene.level))
#创建空白画布
plot(1,xlim=c(1,length(unique(x))+1),ylim=c(1,length(unique(y))+1),
     xaxs="i", yaxs="i",xaxt="n",yaxt="n",
     type="n",bty="n",xlab="",ylab="",
     main = "Coexpression across cancer types",cex.main=2)
#填充颜色
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
#基因名和癌症名
axis(1,at = sort(unique(x)) + 0.5,labels = cancer.level,lty = 0,las = 2)  #添加x轴坐标
axis(2,at = sort(unique(y)) + 0.5,labels = gene.level,lty = 0,las = 1)    #添加y轴坐标
mtext("Cancer types",side = 1,line = 3.5,cex=1.2)    #x轴名称

#绘制图例
par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
barplot(rep(1,length(col.vec2)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec2)),horiz=TRUE,
        axes = F, col=col.vec2)  # Loss
axis(4,at=c(1,ceiling(length(col.vec2)/2),length(col.vec2)),c(round(min(up),1),'Pvalue',round(max(up),1)),tick=FALSE)
par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
barplot(rep(1,length(col.vec1)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec1)),horiz=TRUE,
	axes = F, col=col.vec1)  # Gain
axis(4,at=c(1,ceiling(length(col.vec1)/2),length(col.vec1)),c(round(min(dn),1),'Cor',round(max(dn),1)),tick=FALSE)
dev.off()


