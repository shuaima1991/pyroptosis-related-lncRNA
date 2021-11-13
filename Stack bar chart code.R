setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/46lncRNA census cluster/2分类柱状图/级别柱状图")
library(ggplot2)
#data<-data.frame(Sample<-c(rep('control1',3),rep('control2',3),rep('control3',3),rep('treat1',3),rep('treat2',3),rep('treat3',3),rep('treat4',3)), contion<-rep(c('Cell','Tissue','Organ'),7), value<-c(503,264,148,299,268,98,363,289,208,108,424,353,1,495,168,152,367,146,48,596,143))
data=read.table("输入文件4.txt",sep = "\t",header=T)
#colnames(data)=c('sample',"contion","value")
#Sample是y轴  value是count   fill是图注  x轴填Sample  y轴是
ggplot(data,mapping = aes(Types,Count,fill=Immunetypes))+geom_bar(stat='identity',position='fill') +labs(x = 'Types',y = 'Count') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'))+theme(axis.text.x = element_text(angle = 45, hjust = 1))





