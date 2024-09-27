rm(list=ls())
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(SingleR)
library(celldex)
library(singleseqgset)
library(devtools)
library(grid)
library(gridExtra)
getwd()
dir.create("4-plot")
setwd('4-plot/')
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
sce.all
#Idents(sce.all)

#EPCAM,NKG7,LYZ,CD79A,CLDN5,DCN
marker <- c('EPCAM','NKG7','LYZ','CD79A','CLDN5','DCN')

gene = marker

FeaturePlot(sce.all,features = marker,cols = c("lightgrey" ,"#DE1F1F"),ncol=3,raster=FALSE)
ggsave('FeaturePlot_marker.pdf',width = 12,height = 8)


###使每一个基因的颜色阈值范围调整一致[0-6]
p1 <- FeaturePlot(sce.all, features = marker, combine = FALSE,ncol=3,raster=FALSE )
#colours = c('lightgrey', "#DE1F1F")
fix.sc <- scale_color_gradientn( colours = c('lightgrey', "#DE1F1F"),  limits = c(0, 6))
#+NoLegend()+NoAxes()
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)


##批量画基因,注意图例的范围不同
FeaturePlot(sce.all, features =marker, 
            cols = c("lightgrey", 'red'),
            ncol = 3 ) & NoLegend() & NoAxes() & theme(
              panel.border = element_rect(color = "black", size = 1)
            )

#Featureplot还可以把两个基因画在同一个图中,看右上角可以发现黄色越深的地方两个基因叠加越多
FeaturePlot(sce.all, features = c('S100A9','S100A8'),
            cols = c("lightgrey", "green", "orange"),
            blend=T,blend.threshold=0)

#Featureplot还可以把三个基因画在同一个图中
# 提取tsne坐标
tsne_df <- as.data.frame(sce.all@reductions$umap@cell.embeddings)
tsne_df$cluster <- as.factor(sce.all$celltype)
head(tsne_df)
# 提取基因表达数据并与tsne坐标合并
gene_df <- as.data.frame(GetAssayData(object = sce.all, slot = "data")[c('S100A9','S100A8','CXCL8'), ])

library(ggnewscale)
merged_df <- merge(t(gene_df), tsne_df, by = 0, all = TRUE)
head(merged_df)
colnames(merged_df)
ggplot(merged_df, vars = c("umap_1", "umap_2", 'S100A9','S100A8','CXCL8'), aes(x = umap_1, y = umap_2, colour = S100A9)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("lightgrey", "green"), limits = c(0, 0.3), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = S100A8), size=0.3, alpha=0.7) +
  scale_colour_gradientn(colours = c("lightgrey", "blue"), limits = c(0.1, 0.2), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = CXCL8), size=0.3, alpha=0.1) +
  scale_colour_gradientn(colours = c("lightgrey", "red"), limits = c(0, 0.3), oob = scales::squish)+
  theme_classic()



#2.热图
Idents(sce.all)
table(sce.all$celltype)
Idents(sce.all) = sce.all$celltype
sce1 = sce.all[, sce.all$celltype %in% c( 'B', 'Endothelial','Epithelial', 'Fibro','Myeloid' ,'T&NK' )]

if (!file.exists('sce.markers.csv')) {
  sce.markers <- FindAllMarkers(object = sce1, only.pos = TRUE, 
                                min.pct = 0.25, 
                                thresh.use = 0.25)
  write.csv(sce.markers,file='sce.markers.csv')
} else {
  
  sce.markers = read.csv('sce.markers.csv',row.names = 1)
}


library(dplyr) 
top5 <- sce.markers%>% group_by(cluster) %>% top_n(5, avg_log2FC)

#为了防止数据量太大不好出图，这里在每个亚群提取出来100个
sce.Scale <- ScaleData(subset(sce1,downsample=100),
                       features = top5$gene )
DoHeatmap(sce.Scale,
          features = top5$gene ,
          # group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave('markers_heatmap.pdf',width = 10,height = 7)

top5_dotplot <- DotPlot(sce.all, features = top5$gene)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
top5_dotplot
ggsave('markers_top5_dotplot.pdf',width = 10,height = 7)

setwd('../')
