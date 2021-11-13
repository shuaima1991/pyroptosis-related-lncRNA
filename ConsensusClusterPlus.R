

setwd("")

options(scipen = 1);
rt=read.table("46lncRNA input.txt",header=T,sep="\t",row.names = 1)
d=as.matrix(rt)


colors <- colorRampPalette(c("white", "#FF8C00"))(5)
library(ConsensusClusterPlus)
title=tempdir()
class(d)
results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
                               clusterAlg="hc",distance="pearson",plot="png",tmyPal =colors )




#results[[2]] is theresults result of k=2
results[[4]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]

results[[3]][["consensusClass"]]
c=write.table(results[[3]][["consensusClass"]],"3.txt",sep="\t")
icl = calcICL(results,title=title,plot="png")

icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]

c=results[[2]]$consensusClass
class(c)
c=as.data.frame(c)
b=write.table(icl,"2.txt",sep="\t")
