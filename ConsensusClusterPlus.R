

setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/46lncRNA census cluster")
#取消R科学计数法的代码
options(scipen = 1);
rt=read.table("46lncRNA输入文件.txt",header=T,sep="\t",row.names = 1)
d=as.matrix(rt)

#对上面这个芯片表达数据我们一般会简单的进行normalization （本次采用中位数中心化），然后取在各个样品差异很大的那些gene或者探针的数据来进行聚类分析

colors <- colorRampPalette(c("white", "#FF8C00"))(5)
library(ConsensusClusterPlus)
title=tempdir()
class(d)
results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                               clusterAlg="hc",distance="pearson",plot="png",tmyPal =colors )

#聚类数目K=2，3，4，·····6，采用重抽样方案对样本的80%抽样，经过多次采样，找到稳定可靠的亚组分类。

#然后利用这些有类标签的样本来寻找可以将样本分类的标签基因。可以利用PAM方法找寻标签基因。

#results[[2]] is theresults result of k=2
results[[4]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]
#查k=“”，的样本分类
results[[3]][["consensusClass"]]
c=write.table(results[[3]][["consensusClass"]],"3.txt",sep="\t")
icl = calcICL(results,title=title,plot="png")

icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]
#导出result得list
c=results[[2]]$consensusClass
class(c)
c=as.data.frame(c)
b=write.table(icl,"2.txt",sep="\t")
