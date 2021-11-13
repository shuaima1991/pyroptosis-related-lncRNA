setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/单因素和多因素cox") #设置工作目录

library(forestplot)
#得注意data[,1:3]是不是选择得HR那三列，rt读取文件时rowname=1会错一列一定要注意
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  data=as.matrix(rt)
  HR=data[,1:3]#得注意data[,1:3]是不是选择得HR那三列，rt读取文件时rowname=1会错一列一定要注意
  hr=sprintf("%.3f",HR[,"HR"])
  hrLow=sprintf("%.3f",HR[,"HR.95L"])
  hrHigh=sprintf("%.3f",HR[,"HR.95H"])
  pVal=data[,"pvalue"]
  pVal=ifelse(pVal<0.05, "<0.05", sprintf("%.3f", pVal))
  clrs <- fpColors(box=forestCol,line="darkblue", summary="royalblue")      #定义颜色
  tabletext <- 
    list(c(NA, rownames(HR)),
         append("pvalue", pVal),
         append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   #定义图片文字
  pdf(file=forestFile,width = 9,height = 10,onefile = FALSE)
  forestplot(tabletext, 
             rbind(rep(NA, 3), HR),
             col=clrs,
             graphwidth=unit(50, "mm"),
             xlog=T,
             lwd.ci=4,
             boxsize=0.6,
             xlab="Hazard ratio",
             txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))
  )
  dev.off()
}
############绘制森林图函数############

bioForest(coxFile="multiCox2.txt",forestFile="forest3.pdf",forestCol="red")

