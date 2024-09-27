####14241个上皮细胞，跟文中的12879数量差不多，注释应该大差不差
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
source('scRNA_scripts/lib.R')
source('scRNA_scripts/mycolors.R')
dir.create("7-epi")
setwd("7-epi/") 
getwd()
set.seed(12345)
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
table(sce.all$celltype)
sce1 = sce.all[, sce.all$celltype %in% c( 'Epithelial' )]
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
sce <- FindClusters(sce, resolution = 0.1)
table(sce@meta.data$RNA_snn_res.0.1)  

set.seed(321)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)

sce <- RunUMAP(object = sce, dims = 1:5, do.fast = TRUE)

p = DimPlot(sce,label=T,cols = mycolors) 
p

###先根据文中注释看看情况，再结合CNV注释
genes_to_check = c( #上皮
                   'TM4SF1' ,  
                   'CRABP2',  
                   'UBE2C' , 
                   'TOP2A','MKI67',
                   'CAV1','CLDN18', #AT2
                   'CAPS',  
                   'SCGB1A1') 

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check

p = DotPlot(sce, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p
ggsave('check_markers.pdf',height = 9,width = 7)

##左下角坐标轴
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(sce)

axes <- result$axes
label <- result$label

umap =DimPlot(sce, reduction = "umap",cols = my36colors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.1",label = T,label.box = T) +
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
ggsave('RNA_snn_res.0.1_umap.pdf',width = 9,height = 7)


#####小提琴图
##标记物可视化
VlnPlot(sce,
        features = c('TM4SF1' ,  
                     'CRABP2',  
                     'UBE2C' ,  
                     'CAV1','CLDN18', #AT2
                     'CAPS',  
                     'SCGB1A1'),
        pt.size = 0,
        ncol = 4,
        cols=mycolors)


###如果跟文中对应一下，初步猜测1,3,6,7是正常上皮，其他不确定，有的基因表达和文中对不上

###CNV分析确证####
#1.准备输入文件开始分析之前，需要准备三个文件，https://github.com/broadinstitute/inferCNV/wiki/File-Definitions
#1.1表达矩阵文件。counts矩阵，每列为一个细胞，每行为一个基因。
library(Seurat)
library(infercnv)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
dir.create("CNV")
setwd("CNV") 
epi_sce = sce
Idents(epi_sce)
table(Idents(epi_sce))
###吃亏上当一百年，如果不加前缀，后面0-20直接变成了1-21
epi_sce$cnv <- paste0("C", Idents(epi_sce))
table(epi_sce$cnv)
Idents(epi_sce) = epi_sce$cnv

#抽样，无实际意义仅用于该教程
seurat_object = epi_sce

seurat_object <- subset(seurat_object, downsample=200)
table(Idents(seurat_object))
#筛选分析需要的细胞类型
table(sce.all$celltype)
Endothelial = sce.all[, sce.all$celltype %in% c( 'Endothelial' )]
table(Endothelial$celltype)

T = sce.all[, sce.all$celltype %in% c( 'T&NK' )]
T <- subset(T, downsample=500)

endoMat <-as.data.frame(Endothelial[["RNA"]]$counts )  
TMat <- as.data.frame(T[["RNA"]]$counts) 

spike_mat1 = endoMat
spike_mat2 = TMat

save(spike_mat1,spike_mat2,file = 'reference_mat.Rdata') 

#在还没有运行inferCNV之前不知道上皮细胞恶性与否
epiMat <- as.data.frame(seurat_object[["RNA"]]$counts ) 
epiMat <- epiMat [,unique(colnames(epiMat))]
spike_mat1 = spike_mat1 [,unique(colnames(spike_mat1))]
spike_mat2 = spike_mat2 [,unique(colnames(spike_mat2))]

ids = intersect(rownames(epiMat),rownames(spike_mat1))

dim(spike_mat1[ids,])
dim(spike_mat2[ids,])

this_dat=cbind(epiMat[ids,],spike_mat1[ids,],spike_mat2[ids,])
this_dat[1:6,1:3]

phe = seurat_object@meta.data
length(colnames(this_dat))
length(phe$cnv)
# ref-1和ref-2用作refence group。分别代表Endo和t细胞。spike-1和spike-2分别代表Endo和T细胞，用于掺入epi，评判之后的infercnv的效果。
groupinfo=data.frame(v1=colnames(this_dat),
                     v2=c( phe$cnv ,
                           rep('spike-1',3177),
                           rep('ref-1',1500),
                           rep('spike-2',850),
                           rep('ref-2',500)))
head(groupinfo) 
groupFiles='groupFiles.txt'
table(groupinfo$v2)
write.table(groupinfo,file = groupFiles,
            sep = '\t',quote = F,col.names = F,row.names = F)
print(dim(this_dat))

dat = this_dat
print(dim(dat))
library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),]
dim(dat)
head(dat)[,1:6]
colnames(dat)
str(dat)
expFile='expFile.txt'
rownames(dat)
write.table(dat,file = expFile,sep = '\t',quote = F,row.names = T)
library(data.table)

head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

expFile='expFile.txt'
groupFiles='groupFiles.txt'
geneFile='geneFile.txt'

# duplicate 'row.names' are not allowed
library(infercnv)
table(groupinfo$v2)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c('ref-1',
                                                      'ref-2'))  ## 这个取决于自己的分组信息里面的


# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "infercnv_output",  # dir is auto-created for storing outputs
                              cluster_by_groups=T ,   # cluster
                              hclust_method="ward.D2", plot_steps=T)


save(groupinfo,geneInfor,dat,epi_sce,file = 'infercnv.Rdata')

#应用infercnv的结果
infer_CNV_obj<-readRDS('infercnv_output/run.final.infercnv_obj')
expr<-infer_CNV_obj@expr.data
expr[1:4,1:4]
data_cnv<-as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)

meta = epi_sce@meta.data

###如果从CNV图实在分不清，也可以看看CNVscore
#run.final.infercnv_obj对象文件

tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-1`]
tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-2`]
tmp= cbind(tmp1,tmp2)
down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
oneCopy=up-down
oneCopy
a1= down- 2*oneCopy
a2= down- 1*oneCopy
down;up
a3= up +  1*oneCopy
a4= up + 2*oneCopy 
  
cnv_score_table<-infer_CNV_obj@expr.data
cnv_score_table[1:4,1:4]
cnv_score_mat <- as.matrix(cnv_score_table)
  
# Scoring
cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
# Check
table(cnv_score_table[,1])
# Replace with score 
cnv_score_table_pts <- cnv_score_mat
rm(cnv_score_mat)
# 
cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
cnv_score_table_pts[1:4,1:4]
str(  as.data.frame(cnv_score_table_pts[1:4,1:4])) 
cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  
colnames(cell_scores_CNV) <- "cnv_score" 


head(cell_scores_CNV) 
score=cell_scores_CNV

groupinfo$v1
table(groupinfo$v2 )
a = groupinfo

a$v2 = gsub("ref-1|spike-1", "Endo", a$v2)  
a$v2 = gsub("ref-2|spike-2", "T", a$v2)  

table(groupinfo$v2)
table(a$v2)

identical(groupinfo$v1,rownames(score))
table(score$group)
head(score)
score$group = a$v2
colnames(score)
library(ggthemes)
library(ggplot2)
p <- ggplot(score, aes(x = group, y = cnv_score, fill = group)) +
  geom_boxplot() +theme_base()+   
  theme(legend.position = "none") +
  labs(x = "", y = "CNV Score") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p
###推测1，2，3，5，6有可能是正常上皮


setwd('../')
#####细胞生物学命名
celltype=data.frame(ClusterID=0:8,
                    celltype= 0:8) 


# alveolar type I cell (AT1; AGER+)
# alveolar type II cell (AT2; SFTPA1)
# secretory club cell (Club; SCGB1A1+)
# basal airway epithelial cells (Basal; KRT17+)
# ciliated airway epithelial cells (Ciliated; TPPP3+)
# 


# 这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字
celltype[celltype$ClusterID %in% c(0,5 ),2]='Cancer(TM4SF1+)'
celltype[celltype$ClusterID %in% c( 6 ),2]='Clara-like cancer'
celltype[celltype$ClusterID %in% c(1),2]='AT1'
celltype[celltype$ClusterID %in% c( 8 ),2]='Cancer(CRABP2+)'
celltype[celltype$ClusterID %in% c( 3),2]='AT2'
celltype[celltype$ClusterID %in% c( 2 ),2]='Ciliated cells'
celltype[celltype$ClusterID %in% c(4 ),2]='Cancer(UBE2C+)'
celltype[celltype$ClusterID %in% c( 7 ),2]='Clara cells'

table(epi_sce@meta.data$RNA_snn_res.0.1)
table(celltype$celltype)

epi_sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  epi_sce@meta.data[which(epi_sce@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(epi_sce@meta.data$celltype)


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)
source('../scRNA_scripts/Bottom_left_axis.R')
result <- left_axes(epi_sce)
axes <- result$axes
label <- result$label
celltype_umap =DimPlot(epi_sce, reduction = "umap",cols = my36colors,pt.size = 0.8,
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

saveRDS(epi_sce, "epi_sce_celltype.rds")


Idents(epi_sce) = epi_sce$celltype
if (!file.exists('epimarkers.csv')) {
  markers <- FindAllMarkers(object = epi_sce, only.pos = TRUE, 
                            min.pct = 0.25, 
                            thresh.use = 0.25)
  write.csv(markers,file='epimarkers.csv')
} else {
  
  sce.markers = read.csv('epimarkers.csv',row.names = 1)
}


##比例饼图
source('prop_pie.R')
## 函数说明
# seurat_by：为细胞类型，或seurat_clusters（默认的为细胞类型）
# pheno_by：为样本类型（比如tumor adj）或数据来源等
## provide five palettes
## 1. npg
## 2. jama
## 3. jco
## 4. lancet
## 5. nejm
prop_pie(object = epi_sce,seurat_by = "sample", pheno_by = "celltype",palette = "nejm")

object = epi_sce
seurat_by = "sample"
pheno_by = "celltype"
palette = "nejm"
library(ggrepel)
library(ggsci)
library(tidyverse)
library(ggplot2)
plot_data <- object@meta.data %>% 
  dplyr::select({{pheno_by}}, {{seurat_by}}) 

plot_data <- table(plot_data[[seurat_by]], plot_data[[pheno_by]]) %>%
  as.data.frame() %>%
  group_by(Var2) %>%
  mutate(Total = sum(Freq),
         Proportion = round((Freq / Total)*100,2),
         labels = scales::percent(Freq / Total)) %>%
  mutate(text_y = Total - (Freq/2)) %>%
  dplyr::rename(pheno_by = "Var2") %>%
  dplyr::rename(seurat_by = "Var1")

colors <- colorRampPalette((pal_nejm("default")(8)))(length(unique(plot_data$seurat_by))) #nejm


p <- ggplot(plot_data, aes(x = "", y = Proportion, fill = seurat_by)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) + #theta = y" 意味着使用 y 轴上的值作为角度。因此，在执行 coord_polar(theta = "y") 后，y 轴的值将被解释为角度，并且图形将绘制为极坐标系。
  theme_void() +
  facet_wrap(~ pheno_by) +
  labs(fill = seurat_by) +
  xlab("") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)
p

table(plot_data$seurat_by)
plot_data$seurat_by = factor(plot_data$seurat_by ,levels = c('AIS' ,'MIA' ,'IAC' ))
p2 =  ggplot(plot_data, aes(x = seurat_by, y = Proportion,color = seurat_by)) +
  geom_line(aes(group = 1), color = '#64A36C', size = 1.5) +
  geom_point(size = 5) +
  facet_wrap(~ pheno_by) +
  theme_bw()+
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)
p2
p+p2
###看起来CRABP2+亚群在IAC里骤增
ggsave('pie_prob.pdf',width = 12,height = 7)

setwd('../')

