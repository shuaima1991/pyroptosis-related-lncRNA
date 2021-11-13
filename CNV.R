setwd("D:/Rcode/文章思路/细胞焦亡相关lncRNA/2分类突变全景图/CNV/cluster1 GISTIC")
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
library(BSgenome.Hsapiens.UCSC.hg19)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# Create a chromosomes reference objects function
### segment information


dataSub <- read.table("cluster1.txt",header = T,sep="\t")
chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}

# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg19"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)
#str(chrom)
#pdf("ESCA_copy_number_gistic_score.pdf",12,5)
# Import gistic2 results read gistic output file
scores <- read.table("scores.txt", sep="\t",header=T,stringsAsFactors = F,fileEncoding = "utf-8")
head(scores)
unique(scores$Chromosome)
#把染色体名从阿拉伯数字改为“chr1”、“chrX”的形式
scores[scores$Chromosome==23, "Chromosome"] <- "X"
scores[scores$Chromosome==24, "Chromosome"] <- "Y"
chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",scores$Chromosome))]))

# Important step for accurate length to match back to continual chrom loci
scores$Start.geno <- scores$Start + chrom$chromosomes$chr.length.cumsum[chrID]
scores$End.geno <- scores$End + chrom$chromosomes$chr.length.cumsum[chrID]

# Prepare input data for ploting
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$G.score <- scores.amp$G.score * 1
scores.del <- scores[scores$Type=="Del",]
scores.del$G.score <- scores.del$G.score * -1
scores <- rbind.data.frame(scores.amp,scores.del)
# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim <- c(min(scores$G.score) - 0.1, max(scores$G.score) + 0.1)
title <- paste0("TCGA ESCA overall copy number gistic score", " ", "n=", length(dataSub$patient))

plot(scores.amp$Start.geno, scores.amp$G.score,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gistic score", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$G.score, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .8))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(ylim[1]+0.15,ylim[2]-0.15)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))

#画全部样本的percentage/frequency
#pdf("ESCA_copy_number_percentage.pdf",12,15)
scores.amp <- scores[scores$Type=="Amp",]
scores.amp$frequency<- scores.amp$frequency * 100
scores.del <- scores[scores$Type=="Del",]
scores.del$frequency<- scores.del$frequency * -100

# copy number percentage plot
# seg.col = list(gain = "red", outscale.gain = "darkred", loss = "blue", outscale.red = "midnightblue")
ylim<- c(-90,90)
title=paste0("ESCA overall copy number percentage"," ","n=",length(dataSub$patient))

plot(scores.amp$Start.geno, scores.amp$frequency,
     pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i", 
     xlim = c(0,chrom$genome.length), ylim = ylim,
     main = title, cex.main = 2, ylab = "gain/loss percentage in cohort", xlab = NA,
     cex.lab = 2, col = adjustcolor("darkred", alpha.f = .8), xaxt = "n", lwd = 2, las=1) # las=1 rotating axis labels in R
lines(scores.del$Start.geno, scores.del$frequency, type='h', lwd = 2, col = adjustcolor("midnightblue", alpha.f = .5))
ink <- chrom$chromosomes$chrN %in% chrID
yrange = abs(diff(ylim))
m.pos <- c(-80,80)
m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
try(text(x = chrom$chromosomes$middle.chr.geno[ink], y = m.pos[m.mod], labels = chrom$chromosomes$chrom[ink], cex = 1))
abline(h = 0.0, col = 1, lwd = 1, lty = 3)
abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)

col1 <- adjustcolor("darkred", alpha.f = .8)
col2 <- adjustcolor("midnightblue", alpha.f = .8)
# The position of the legend can be specified also using the following keywords : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))

