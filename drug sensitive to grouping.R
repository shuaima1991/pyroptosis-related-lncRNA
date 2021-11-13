options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/TCGA高低分组/药物反应性")
library(pRRophetic)
library(ggplot2)
library(cowplot)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#表达矩阵
dat <- read.table("输入文件.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
dat=dat[!duplicated(dat$id), ]
rownames(dat)=dat[,1]
dat=dat[,-1]
dat[1:3, 1:3]
#分组信息
ann <- read.table("group.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
head(ann)
table(ann$ImmClust)
#药物名字
GCP.drug <- read.table("drug_eg.txt") #如果要例文的两种药物，就换成drug_eg.txt
GCP.drug <- GCP.drug$V1
#这里以前12种药物为例
GCP.drug <- GCP.drug[1:28]
# 自定义足够多的box的颜色，颜色数量至少等于分组数量
jco <- c("#EABF00", "#2874C5", "red")

### 药物敏感性预测 ###、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()

for (drug in GCP.drug) {
  set.seed(1248103) # 因为预测过程默认10-fold CV，所以设置种子以便结果可重复
  cat(drug," starts!\n") # 提示当前药物已开始分析
  
  # 预测IC50值，采用默认参数，详细可参考??pRRopheticPredict参数列表
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) # 1表示若有重复基因取均值处理
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} # 若名字不匹配则报错退出
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "ImmClust"=ifelse(ann$ImmClust == "C1","highrisk","lowrisk"), # 这里我修改了C1和C2的名字
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("highrisk","lowrisk"),ordered = T) # 把类改成因子变量
  
  # 绘图
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=ImmClust, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = ImmClust)) + 
    scale_fill_manual(values = jco[1:length(unique(ann$ImmClust))]) + #自定义box的配色
    theme(legend.position="none") + # 倾斜字体
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) # 补上title
  
  plotp[[drug]] <- p # 保存在列表里供合并图片用
  cat(drug," has been finished!\n") # 提示当前药物已分析结束
}

# 合并图片
#适合展示两种药物
p1 <- plot_grid(plotp[[1]],plotp[[2]],labels = c("A","B"),nrow = 1) # title可以AI下拉到合适位置，就如例文所示
ggsave("boxplot of predicted IC50.pdf", width = 6, height = 5)

# 适合展示多种药物
p2 <- plot_grid(plotlist=plotp, ncol=6)
ggsave("boxplot of predicted IC50_multiple.pdf", width = 15, height = 12)

p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "highrisk"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "lowrisk"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) # 两组样本秩和检验p值
}
names(p) <- GCP.drug
print(p) #打印p值，因为这里有一个不显著，所以当时没有放在boxplot上，有需要的也可以加在ggplot的title里，或参考FigureYa12box直接画在图上。
#保存到文件
write.table(p,"output_pvalue3.txt", quote = F, sep = "\t")

