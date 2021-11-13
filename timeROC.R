setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/TCGA高低分组/时间ROC")
rt=read.table("输入文件.txt",header=T,sep="\t",check.names=F)    #读取cox回归风险文件
library(timeROC)
library(survival)
rt$futime=rt$futime/12
ROC_rt <- timeROC(T=rt$futime,delta=rt$fustat,

                 marker=rt$riskScore,cause=1,

                 weighting='marginal',

                 times=c(1,3,5,7,10),ROC=TRUE)
#画图

plot(ROC_rt,time=1,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='green',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=7,col='yellow',add=TRUE,title=FALSE,lwd=2)

plot(ROC_rt,time=10,col='purple',add=TRUE,title=FALSE,lwd=2)

legend('bottomright',

      c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],2)),
        paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],2)),
        paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],2)),
        paste0('AUC at 7 years: ',round(ROC_rt$AUC[4],2)),

        paste0('AUC at 10 years: ',round(ROC_rt$AUC[5],2))),

      col=c('red','blue',"green",'yellow',"purple"),lwd=2,bty = 'n')

