rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
source('scRNA_scripts/lib.R')
source('scRNA_scripts/mycolors.R')
dir.create("12-T")
setwd("12-T/") 
getwd()
set.seed(12345)
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
table(sce.all$celltype)
sce1 = sce.all[, sce.all$celltype %in% c( 'T&NK' )]
LayerData(sce1, assay = "RNA", layer = "counts")
sce1 <- JoinLayers(sce1)

as.data.frame(sce1@assays$RNA$counts[1:10, 1:2])
head(sce1@meta.data, 10)
table(sce1$orig.ident) 

sce = sce1
sce <- NormalizeData(sce, normalization.method =  "LogNormalize",  
                     scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000)  
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce) 

sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.3)
table(sce@meta.data$RNA_snn_res.0.3)  

set.seed(321)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)

sce <- RunUMAP(object = sce, dims = 1:5, do.fast = TRUE)

p = DimPlot(sce,label=T,cols = mycolors) 
p

###先根据文中注释看看情况
genes_to_check = c( 'CD4',
                   'LTB' ,   #memory CD4
                   'CXCL13',  #exhausted CD4
                   'CD8',
                   'GZMK' ,  # CD8
                   'CCL4L2', # CD8
                   'GNLY', #NK
                   'NKG7','FGFBP2' , #NK
                   'HSPA1B', #Fibro-like
                   'FOXP3') #Treg

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

p = DotPlot(sce, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p
ggsave('check_markers.pdf',height = 10,width = 7)

##左下角坐标轴
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(sce)

axes <- result$axes
label <- result$label

umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.3",label = T,label.box = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
umap
ggsave('RNA_snn_res.0.3_umap.pdf',width = 9,height = 7)



#####细胞生物学命名
celltype=data.frame(ClusterID=0:11,
                    celltype= 0:11) 


# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
celltype[celltype$ClusterID %in% c(0,1 ),2]='memory_CD4'
celltype[celltype$ClusterID %in% c( 10 ),2]='exhausted_CD4'
celltype[celltype$ClusterID %in% c(4),2]='GZMK_CD8'
celltype[celltype$ClusterID %in% c( 8 ),2]='CCL4L2_CD8'
celltype[celltype$ClusterID %in% c(3,6,11 ),2]='GNLY_NK'
celltype[celltype$ClusterID %in% c( 2 ),2]='Fibro-like T'
celltype[celltype$ClusterID %in% c(7,9 ),2]='FGFBP2_NK'
celltype[celltype$ClusterID %in% c( 5),2]='Treg'

table(sce@meta.data$RNA_snn_res.0.3)
table(celltype$celltype)

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)

result <- left_axes(sce)
axes <- result$axes
label <- result$label
celltype_umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.5,
                       group.by = "celltype",label = T) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
celltype_umap
ggsave('umap_by_celltype.pdf',width = 9,height = 7)

saveRDS(sce, "Tsce_celltype.rds")



setwd('../')

