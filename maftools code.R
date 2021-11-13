setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/2分类突变全景图/突变瀑布图")
library(maftools)
luad <- read.maf(maf="TCGA.GBM.LGG合并后.maf")

# 从临床数据中提取性别对应的"Tumor_Sample_Barcode"
clin <- read.table("TCGA-GBM.LGG_phenotype .tsv", header=T, sep="\t")
clin.IS1 <- subset(clin, Type=="cluster1")$Tumor_Sample_Barcode
clin.IS2 <- subset(clin, Type=="cluster2")$Tumor_Sample_Barcode

# 使用subsetMaf构建男性和女性的MAF对象
luad.IS1 <- subsetMaf(maf=luad, tsb=clin.IS1, isTCGA=TRUE)
luad.IS2 <- subsetMaf(maf=luad, tsb=clin.IS2, isTCGA=TRUE)

hnsc.titv = titv(maf = luad.IS1, plot = FALSE, useSyn = TRUE)
#plot titv summary
 plotTiTv(res = hnsc.titv)
hnsc.titv = titv(maf = luad.IS2, plot = FALSE, useSyn = TRUE)
 #plot titv summary
 plotTiTv(res = hnsc.titv)

#drive基因
shnsc.sigIS1 = oncodrive(maf = luad.IS1, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS1, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
shnsc.sigIS2 = oncodrive(maf = luad.IS2, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = shnsc.sigIS2, fdrCutOff = 0.01, useFraction = TRUE,labelSize = 1)
#VAF
 plotVaf(maf = luad.IS1, vafCol = NULL)
 plotVaf(maf = luad.IS2, vafCol = NULL)
 #此处使用RColorBrewer的颜色，当然也可以使用任意颜色
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
#查看变异类型对应的颜色
print(vc_cols)
#>   Frame_Shift_Del Missense_Mutation Nonsense_Mutation         Multi_Hit 
#>         "#A6CEE3"         "#1F78B4"         "#B2DF8A"         "#33A02C" 
#>   Frame_Shift_Ins      In_Frame_Ins       Splice_Site      In_Frame_Del 
#>         "#FB9A99"         "#E31A1C"         "#FDBF6F"         "#FF7F00"

oncoplot(luad.IS1, colors = vc_cols, top = 20)
oncoplot(luad.IS2, colors = vc_cols, top = 20)

output <- somaticInteractions(maf=luad.IS1, top=25, pvalue=c(0.05, 0.01))
output <- somaticInteractions(maf=luad.IS2, top=25, pvalue=c(0.05, 0.01))

output <- somaticInteractions(maf=luad,genes = c("CASP1","CASP3","CASP4","CASP5","CASP6","CASP7","CASP8","GSDMB","GSDMC","GSDMD","GSDME","GZMA","GZMB"), pvalue=c(0.05, 0.01))
#VAF
plotVaf(maf = luad.IS1,width = 100,height = 10)
plotVaf(maf = luad.IS2,width = 100,height = 10)
plotVaf(maf = luad.IS3,width = 100,height = 10)
plotVaf(maf = luad.IS4,width = 100,height = 10)
luad.vaf <- vafCompare(m1 = luad.IS2,m2 = luad.IS3)
#ATCG
hnsc.titv = titv(maf = luad, plot = FALSE, useSyn = TRUE)
write.table(hnsc.titv[["raw.counts"]],"ATCG.txt",sep="\t")
# 使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=luad.IS1, m2=luad.IS2, m3=luad.IS3, m4=luad.IS4, m1Name="IS1", m2Name="IS2",m3Name="IS3",m4Name="IS4", minMut=5)
# 结果保存到文件"female_vs_male.tsv"
write.table(fvsm$results, file="female_vs_male.tsv", quote=FALSE, row.names=FALSE, sep="\t")


#森林图
fvsm1 <- mafCompare(m1=luad.IS1, m2=luad.IS2, m1Name="IS1", m2Name="IS2", minMut=5)
forestPlot(mafCompareRes=fvsm1, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
fvsm2 <- mafCompare(m1=luad.IS3, m2=luad.IS4, m1Name="IS3", m2Name="IS4", minMut=5)
forestPlot(mafCompareRes=fvsm2, pVal=0.001, color=c("maroon", "royalblue"), geneFontSize=0.8)
#survival
mafSurvival(maf=luad.IS1, clinicalData = clin,genes=c("IDH1"), time="days_to_death.demographic", Status="vital_status.demographic", isTCGA=TRUE)
mafSurvival(maf=luad.IS2, clinicalData = clin,genes=c("IDH1"), time="days_to_death.demographic", Status="vital_status.demographic", isTCGA=TRUE)


lollipopPlot(maf = luad.IS1, gene = 'IDH1', AACol = 'HGVSp_Short', showMutationRate = TRUE)



lollipopPlot(m1 = luad.IS1, m2 = luad.IS2, gene = "IDH1", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "luad.IS1", m2_name = "luad.IS2")

#导出文件
write.table(luad@variants.per.sample, file="TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad@variant.type.summary, file="variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad@variant.classification.summary, file="variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

write.table(luad.IS2@variants.per.sample, file="IS2 TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS2@variant.type.summary, file="IS2 variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS2@variant.classification.summary, file="IS2 variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

write.table(luad.IS3@variants.per.sample, file="IS3 TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS3@variant.type.summary, file="IS3 variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS3@variant.classification.summary, file="IS3 variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

write.table(luad.IS4@variants.per.sample, file="IS4 TMB.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS4@variant.type.summary, file="IS4 variant.type.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(luad.IS4@variant.classification.summary, file="IS4 variant.classification.txt", quote=FALSE, row.names=FALSE, sep="\t")

