
pFilter=0.05                                                      #定义单因素显著性

setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/单因素和多因素cox") #设置工作目录
library(survival)                                                 #引用包
library(UpSetR)
rt=read.table("单因素输入文件.txt",header=T,sep="\t",check.names=F,row.names=1)       #读取输入文件

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 outTab=rbind(outTab,
              cbind(id=i,
                    z=coxSummary$coefficients[,"z"],
                    HR=coxSummary$conf.int[,"exp(coef)"],
                    HR.95L=coxSummary$conf.int[,"lower .95"],
                    HR.95H=coxSummary$conf.int[,"upper .95"],
                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
              )
}
#输出所有单因素的结果
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著的结果
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著AS的PSI值，用于后续建模
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


#绘制upset图
gene=sapply(strsplit(sigGenes,"\\|"),"[",1)
asType=sapply(strsplit(sigGenes,"\\|"),"[",3)
upsetList=list(AA=unique(gene[asType=="AA"]),
               AD=unique(gene[asType=="AD"]),
               AP=unique(gene[asType=="AP"]),
               AT=unique(gene[asType=="AT"]),
               ES=unique(gene[asType=="ES"]),
               ME=unique(gene[asType=="ME"]),
               RI=unique(gene[asType=="RI"]) )
upsetData=fromList(upsetList)

pdf(file="uniCoxUpset.pdf",onefile = FALSE,width=8,height=5)              #保存图片
upset(upsetData,
      nsets = 7,                                    #展示可变剪切类型个数
      order.by = "freq",                            #按照数目排序
      show.numbers = "yes",                         #柱状图上方是否显示数值
      number.angles = 20,                           #字体角度
      point.size = 1.5,                             #点的大小
      matrix.color="red",                           #交集点颜色
      line.size = 0.8,                              #线条粗线
      mainbar.y.label = "Gene Intersections",
      sets.x.label = "Set Size")
dev.off()

######Video source: http://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056