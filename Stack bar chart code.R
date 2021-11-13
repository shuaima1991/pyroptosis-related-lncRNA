setwd("")
library(ggplot2)
#data<-data.frame(Sample<-c(rep('control1',3),rep('control2',3),rep('control3',3),rep('treat1',3),rep('treat2',3),rep('treat3',3),rep('treat4',3)), contion<-rep(c('Cell','Tissue','Organ'),7), value<-c(503,264,148,299,268,98,363,289,208,108,424,353,1,495,168,152,367,146,48,596,143))
data=read.table("input.txt",sep = "\t",header=T)
#colnames(data)=c('sample',"contion","value")

ggplot(data,mapping = aes(Types,Count,fill=Immunetypes))+geom_bar(stat='identity',position='fill') +labs(x = 'Types',y = 'Count') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'))+theme(axis.text.x = element_text(angle = 45, hjust = 1))





