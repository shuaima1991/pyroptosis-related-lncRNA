setwd("")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
rt=read.table("46lncRNAinput.txt",header = T,row.names = 1,sep = "\t")
metadata=read.table("metadata.tsv.txt",header = T,row.names = 1,sep="\t")
rt=as.matrix(rt)

metadata[1:6,1:4]
Type <- metadata$Type
Type <- Type[!is.na(Type)] 
length(unique(Type)) #5
# 
pick.col <- brewer.pal(9, 'Greens') # allowed maximum for palette Greens is 9
col.Type <- colorRampPalette(pick.col)(length(unique(Type)))

# Grade
Grade <- metadata$Grade
Grade <- Grade[!is.na(Grade)]
pick.col <- brewer.pal(9, 'Blues')
col.Grade <- colorRampPalette(pick.col)(length(unique(Grade)))

# Gender
Gender <- metadata$Gender
Gender <- Gender[!is.na(Gender)]
pick.col <- brewer.pal(9, 'Oranges')
col.Gender <- colorRampPalette(pick.col)(length(unique(Gender)))

# Original.Subtype
Original.Subtype <- metadata$Original.Subtype
Original.Subtype <- Original.Subtype[!is.na(Original.Subtype)]
pick.col <- brewer.pal(9, 'Purples')
col.Original.Subtype <- colorRampPalette(pick.col)(length(unique(Original.Subtype)))
#
IDH.status <- metadata$IDH.status
IDH.status <- IDH.status[!is.na(IDH.status)]
pick.col <- brewer.pal(9, 'Purples')
col.IDH.status <- colorRampPalette(pick.col)(length(unique(IDH.status)))
#
Vital.status <- metadata$Vital.status
Vital.status <- Vital.status[!is.na(Vital.status)]
pick.col <- brewer.pal(9, 'Purples')
col.Vital.status <- colorRampPalette(pick.col)(length(unique(Vital.status)))
#
immune.score<- metadata$immune.score
immune.score <- immune.score[!is.na(immune.score)]
pick.col <- brewer.pal(9, 'Purples')
col.immune.score <- colorRampPalette(pick.col)(length(unique(immune.score)))
#
stromal.score<- metadata$stromal.score
stromal.score <- stromal.score[!is.na(stromal.score)]
pick.col <- brewer.pal(9, 'Purples')
col.stromal.score <- colorRampPalette(pick.col)(length(unique(stromal.score)))

# 
ann <- data.frame(
  Type = metadata$Type,
  Grade = metadata$Grade,
  Gender = metadata$Gender,
  Original.Subtype = metadata$Original.Subtype,
  IDH.status = metadata$IDH.status,
  Vital.status=metadata$Vital.status,
  immune.score= metadata$immune.score,
  stromal.score= metadata$stromal.score)
# 
names(col.Type)=as.character(0:1)
names(col.Grade)=as.character(0:2)
names(col.Gender)=as.character(0:1)
names(col.Original.Subtype)=as.character(0:7)
names(col.IDH.status)=as.character(0:1)
names(col.Vital.status)=as.character(0:1)
names(col.immune.score)=as.character(0:4)
names(col.stromal.score)=as.character(0:4)
colors=list(Type = c('cluster1' = 'blue', 'cluster2' = 'red'),
            Grade=c("G4"="yellow","G3"="brown","G2"="purple"),
            Gender=c("female"="orange","male"="tan"),
            Original.Subtype=c("Proneural"="violet","Mesenchymal"="brown","G-CIMP"="blue","Classical"="green","Neural"="red","IDHwt"="magenta","IDHmut-non-codel"="navy","IDHmut-codel"="purple"),
            IDH.status=c("Mutant"="green","WT"="red"),
            Vital.status=c("1"="turquoise", "0"="coral"),
            immune.score=c("1"="#FFCCCC","2"="#FF8888","3"="#FF3333","4"="#FF0000","5"="#CC0000"),
            stromal.score=c("1"="#33FFFF","2"="#5555FF","3"="#9955FF","4"="#E93EFF","5"="#FF3EFF"))
colAnn <- HeatmapAnnotation(
  df = ann,
  col = colors,
  which = 'col', # set 'col' (samples) or 'row' (gene) annotation
  na_col = 'white', # NA
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  # 
  annotation_legend_param = list(
    Type = list(
      nrow = 2,
      title = 'Type',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Grade = list(
      nrow = 3,
      title = 'Grade',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Gender = list(
      nrow = 2,
      title = 'Gender',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Original.Subtype = list(
      nrow = 8,
      title = 'Original.Subtype',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    IDH.status = list(
      nrow = 2,
      title = 'IDH.status',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    immune.score = list(
      nrow = 7,
      title = 'immune.score',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    stromal.score = list(
      nrow = 5,
      title = 'stromal.score',
      title_position = 'topcenter',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold'))))
# 
library(circlize)
rt=t(scale(t(rt)))
col_fun = colorRamp2(c(-1, 0, 1),c('#1E90FF', 'white', '#FF4500'))
group_list=c(rep('cluster1',412),rep('cluster2',197))
#cluster_columns = FALSE
Heatmap(rt, name = "mat", col = col_fun,show_column_names = FALSE,top_annotation = colAnn,cluster_columns = FALSE,column_split = group_list)

