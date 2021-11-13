
setwd("") 

library(survival)
library(survminer)

rt=read.table("input.txt",header=T,sep="\t",check.names=F,row.names=1)
#rt[,"VCAN"]=log2(rt[,"VCAN"]+1)

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
    width = 7,            
    height = 6,            
)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()
