
setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/单因素和多因素cox") #设置工作目录

library(survival)
library(survminer)

rt=read.table("多因素输入文件.txt",header=T,sep="\t",check.names=F,row.names=1)
#rt[,"VCAN"]=log2(rt[,"VCAN"]+1)
#.代表其他所有列
rt$Survival_month=rt$Survival_month+0.04
multiCox=coxph(Surv(rt$Survival_month, rt$status) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox2.xls",sep="\t",row.names=F,quote=F)

pdf(file="forest.pdf",
    width = 7,             #图片的宽度
    height = 6,            #图片的高度
)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()
