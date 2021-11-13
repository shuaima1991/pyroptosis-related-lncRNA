
library("glmnet")
library("survival")

setwd("")               
rt=read.table("lasso输入文件.txt",header=T,sep="\t",row.names=1)       
#rt$Survival_time=rt$Survival_time/365
rt$Survival_time=rt$Survival_time+0.04
rt$Survival_time=rt$Survival_time/12

gene=read.table("gene.txt",header=F)
#rt=rt[,c("Survival_time","Vital_status",as.vector(gene[,1]))]

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$Survival_time,rt$Vital_status))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
outCol=c("Survival_time","Vital_status",lassoGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(riskScore),risk)
write.table(cbind(id=rownames(outTab),outTab),
    file="lassoRisk.txt",
    sep="\t",
    quote=F,
    row.names=F)
