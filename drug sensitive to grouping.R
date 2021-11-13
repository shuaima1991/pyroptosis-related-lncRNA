options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
setwd(")
library(pRRophetic)
library(ggplot2)
library(cowplot)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
dat <- read.table("input.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
dat=dat[!duplicated(dat$id), ]
rownames(dat)=dat[,1]
dat=dat[,-1]
dat[1:3, 1:3]
#分组信息
ann <- read.table("group.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
head(ann)
table(ann$ImmClust)

GCP.drug <- read.table("drug_eg.txt") #如果要例文的两种药物，就换成drug_eg.txt
GCP.drug <- GCP.drug$V1

GCP.drug <- GCP.drug[1:28]

jco <- c("#EABF00", "#2874C5", "red")


GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()

for (drug in GCP.drug) {
  set.seed(1248103) 
  cat(drug," starts!\n") 
  
  
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) 
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")}
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "ImmClust"=ifelse(ann$ImmClust == "C1","highrisk","lowrisk"), 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("highrisk","lowrisk"),ordered = T) 
  

  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=ImmClust, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = ImmClust)) + 
    scale_fill_manual(values = jco[1:length(unique(ann$ImmClust))]) + 
    theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) 
  
  plotp[[drug]] <- p 
  cat(drug," has been finished!\n") 
}

p1 <- plot_grid(plotp[[1]],plotp[[2]],labels = c("A","B"),nrow = 1) # title可以AI下拉到合适位置，就如例文所示
ggsave("boxplot of predicted IC50.pdf", width = 6, height = 5)


p2 <- plot_grid(plotlist=plotp, ncol=6)
ggsave("boxplot of predicted IC50_multiple.pdf", width = 15, height = 12)

p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "highrisk"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$ImmClust %in% "lowrisk"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp) 
}
names(p) <- GCP.drug
print(p) 
write.table(p,"output_pvalue3.txt", quote = F, sep = "\t")

