###### step2: 确定单细胞亚群生物学名字 ######
# 一般来说，为了节省工作量，我们选择0.1的分辨率进行命名
# 因为命名这个步骤是纯人工操作
# 除非0.1确实分群太粗狂了，我们就选择0.8 

###### 常见分群
# T Cells (CD3D, CD3E, CD8A), 
# B cells (CD19, CD79A, MS4A1 [CD20]), 
# Plasma cells (IGHG1, MZB1, SDC1, CD79A), 
# macrophages (CD68, CD163),
# 'CCL3L1' ,   #M2
# 'FABP4',  #M1
# Monocytes  (CD14),
# NK Cells (FGFBP2, FCG3RA, CX3CR1),  
# Photoreceptor cells (RCVRN), 
# Fibroblasts (FGF7, MME), 
# Neutrophil ('G0S2', 'S100A9','S100A8','CXCL8')
# Endothelial cells (PECAM1, VWF). 
# epi or tumor (EPCAM, KRT19, PROM1, ALDH1A1, CD24).
# immune (CD45+,PTPRC), epithelial/cancer (EpCAM+,EPCAM), 
# stromal (CD10+,MME,fibo or CD31+,PECAM1,endo) 


#####注释亚群
#通常我们第一层次降维聚类分群：immune (CD45+,PTPRC),epithelial/cancer (EpCAM+,EPCAM),stromal (CD10+,MME,fibro or CD31+,PECAM1,endo)
#文中提供的都是常见的细胞群：上皮细胞（EPCAM、KRT19、CLDN4）、基质（PECAM1、CLO1A2、VWF）、增殖性（MKI67、STMN1、PCNA）、T（CD3D、CD3E、CD2）、B（CD79A，IGHG1，MS4A1），NK（KLRD1、GNLY、KLRF1）和髓系（CSF1R、CSF3R、CD68）细胞。
rm(list=ls())
source('scRNA_scripts/lib.R')
source('scRNA_scripts/mycolors.R')
sce.all.int = readRDS('2-harmony/sce.all_int.rds') 
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
colnames(sce.all.int@meta.data) 

dir.create("./3-Celltype")
setwd("./3-Celltype")
scRNA=sce.all.int

genes_to_check = c('EPCAM','KRT19','CLDN4','SCGB1A1',  #上皮
                   'PECAM1' , 'CLO1A2', 'VWF',  #基质
                   'CDH5', 'PECAM1', 'VWF','CLDN5',  #内皮
                   'LUM' , 'FGF7', 'MME',  #成纤维
                   'CD3D', 'CD3E', 'CD8A', 'CD4','CD2', #T
                   'AIF1', 'C1QC','C1QB','LYZ',  #巨噬
                   'MKI67', 'STMN1', 'PCNA',  #增殖
                   'CPA3' ,'CST3', 'KIT', 'TPSAB1','TPSB2',#肥大
                   'GOS2', 'S100A9','S100A8','CXCL8', #中性粒细胞
                   'KLRD1', 'GNLY', 'KLRF1','AREG', 'XCL2','HSPA6', #NK
                   'MS4A1','CD19', 'CD79A','IGHG1','MZB1', 'SDC1',  #B
                   'IGHD',  #MALT B
                   'CSF1R', 'CSF3R', 'CD68') #髓系

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

p = DotPlot(scRNA, features = unique(genes_to_check),
                assay='RNA'  )  + coord_flip()
p
#ggsave('check_last_markers.pdf',height = 11,width = 11)

####构建左下角坐标轴
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(scRNA)
axes <- result$axes
label <- result$label


umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
                  group.by = "RNA_snn_res.0.5",label = T,label.box = T) +
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
ggsave('RNA_snn_res.0.5_umap.pdf',width = 9,height = 7)

##图中的标签和方框都是可以自定义的，例如下面这副删掉label和label.box
sample_umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
                     group.by = "sample") +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
  theme(plot.title = element_blank())
sample_umap
ggsave('sample_umap.pdf',width = 9,height = 7)

umap+sample_umap


#####细胞生物学命名
celltype=data.frame(ClusterID=0:22,
                    celltype= 0:22) 

# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
celltype[celltype$ClusterID %in% c(3,6,9,13,14,15,16,20,22 ),2]='Myeloid'
celltype[celltype$ClusterID %in% c( 5,11,12,17,19 ),2]='Epithelial'
celltype[celltype$ClusterID %in% c(0,1,21),2]='T&NK'
celltype[celltype$ClusterID %in% c( 8 ),2]='Fibro'
celltype[celltype$ClusterID %in% c( 18 ),2]='Proliferative'
celltype[celltype$ClusterID %in% c( 10 ),2]='Plasma'
celltype[celltype$ClusterID %in% c( 2),2]='B'
celltype[celltype$ClusterID %in% c( 7 ),2]='Endothelial'
celltype[celltype$ClusterID %in% c( 4),2]='Mast'

table(scRNA@meta.data$RNA_snn_res.0.5)
table(celltype$celltype)

scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)
celltype_umap =DimPlot(scRNA, reduction = "umap",cols = my36colors,pt.size = 0.8,
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

umap + sample_umap+celltype_umap 
ggsave('combine_umap.pdf',width = 15,height = 7)

saveRDS(scRNA, "sce_celltype.rds")

setwd('../')



