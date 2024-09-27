rm(list=ls())
source('scRNA_scripts/lib.R')
dir.create("6-DEG")
setwd('6-DEG/')
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
sce.all
library(tidyverse)
library(tinyarray)
scRNA = sce.all
head(scRNA@meta.data)
Idents(scRNA) = scRNA$sample
ct <- c('IAC', 'MIA', 'AIS')


if (!file.exists('markers_list.Rdata')) {
  markers_list <- list()
  for (i in 1:(length(ct)-1)) {
    for (j in (i+1):length(ct)) {
      markers <- FindMarkers(scRNA, group.by = "sample",
                             logfc.threshold = 0.1,
                             ident.1 = ct[i],
                             ident.2 = ct[j])
      markers_list[[paste0(ct[i], "_vs_", ct[j])]] <- markers
    }
  }
  save(markers_list,file = 'markers_list.Rdata')
} else {
  
  load('markers_list.Rdata')
}

length(markers_list)
lapply(markers_list,nrow)

all_markers_sig = lapply(markers_list, function(x){
  markers_sig <- subset(x, p_val_adj < 0.01)
})


marker_stat = as.data.frame(lapply(all_markers_sig,function(x){
  # x=all_markers_sig[[1]]
  Up = sum(x$avg_log2FC>1)
  Down = sum(x$avg_log2FC< -1)
  Total = Up+Down
  return(c(Up, Down, Total))
}))
rownames(marker_stat) = c("Up","Down","Total")
marker_stat
library(gridExtra)

IAC_vs_MIA = all_markers_sig[[1]]
IAC_vs_AIS = all_markers_sig[[2]]
MIA_vs_AIS = all_markers_sig[[3]]
##每一组加上group上下调
deg <- function(dat){
  dat$group = NA
  dat$group[(dat$avg_log2FC < -1) & (dat$p_val_adj < 0.01)] <- "down"
  dat$group[(dat$avg_log2FC > 1) & (dat$p_val_adj < 0.01)] <- "up"
  table(dat$group )
  dat = na.omit(dat)
  dat
}
P001 = deg(IAC_vs_MIA)
table(P001$group)
P002 = deg(IAC_vs_AIS)
P003 = deg(MIA_vs_AIS)

####veen图####
###上调
IAC_vs_MIA_up = rownames(P001[P001$group == 'up',]) 
IAC_vs_AIS_up = rownames(P002[P002$group == 'up',]) 
MIA_vs_AIS_up = rownames(P003[P003$group == 'up',]) 

inter_upgene <- intersect(intersect(IAC_vs_MIA_up, IAC_vs_AIS_up), MIA_vs_AIS_up)


#三元#
#BiocManager::install("VennDetail")
library(VennDetail)
library(VennDiagram) 
#IAC_vs_MIA_up,IAC_vs_AIS_up,MIA_vs_AIS_up
venn <- venndetail(list(IAC_vs_MIA_up = IAC_vs_MIA_up, IAC_vs_AIS_up = IAC_vs_AIS_up, MIA_vs_AIS_up= MIA_vs_AIS_up))
detail(venn) 

# 韦恩图
venn.diagram(x=list(IAC_vs_MIA_up,IAC_vs_AIS_up,MIA_vs_AIS_up),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             category.names = c("IAC_vs_MIA_up", "IAC_vs_AIS_up","MIA_vs_AIS_up") , #标签名
             cat.dist = 0.07, # 标签距离圆圈的远近
             cat.pos = c(-30, -330, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 1, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='veenup.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw",# 压缩算法
             height = 2100 ,   # 高度
             width = 2300
             
)
###结果和文章中不一致，没关系继续往下

# 韦恩饼图
plot(venn, type = "vennpie")

##多个数据集的韦恩饼图
# vennpie(venn, 
#         min = 4 # 显示集合至少包含来自四个数据集的元素
#         # any = 1, revcolor = "lightgrey" # 突出显示唯一或共享子集
# )

# 韦恩条形图
dplot(venn, order = TRUE, textsize = 4)

# upset图
plot(venn, type = "upset")


###下调
IAC_vs_MIA_down = rownames(P001[P001$group == 'down',]) 
IAC_vs_AIS_down = rownames(P002[P002$group == 'down',]) 
MIA_vs_AIS_down = rownames(P003[P003$group == 'down',]) 

inter_downgene <- intersect(intersect(IAC_vs_MIA_down, IAC_vs_AIS_down), MIA_vs_AIS_down)

venn <- venndetail(list(IAC_vs_MIA_down = IAC_vs_MIA_down, IAC_vs_AIS_down = IAC_vs_AIS_down, MIA_vs_AIS_down= MIA_vs_AIS_down))
detail(venn) 

# 韦恩图
venn.diagram(x=list(IAC_vs_MIA_down,IAC_vs_AIS_down,MIA_vs_AIS_down),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             category.names = c("IAC_vs_MIA_down", "IAC_vs_AIS_down","MIA_vs_AIS_down") , #标签名
             cat.dist = 0.07, # 标签距离圆圈的远近
             cat.pos = c(-30, -330, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex = 1, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='veendown.png',# 文件保存
             imagetype="png",  # 类型（tiff png svg）
             resolution = 400,  # 分辨率
             compression = "lzw",# 压缩算法
             height = 2100 ,   # 高度
             width = 2300
             
)
###结果和文章中不一致，没关系继续往下

# 韦恩饼图
plot(venn, type = "vennpie")

# 韦恩条形图
dplot(venn, order = TRUE, textsize = 4)

# downset图
plot(venn, type = "upset")

##3D PCA
set.seed(98765)
##PCA showing the DEGs among the three groups. The distance between dots represents the difference between groups
table(scRNA$patient)
meta = scRNA@meta.data
meta$sample_patient = paste(meta$sample,meta$patient,sep = '_')
table(meta$sample_patient)
#Idents(scRNA)
sce = scRNA
colnames(sce)
rownames(sce)
sce@meta.data = meta
avg <- AverageExpression(object = sce, group.by = "sample_patient")
a = avg$RNA
b = as.data.frame(a)
###提取出差异基因
DEGgene = c(rownames(IAC_vs_AIS),rownames(IAC_vs_MIA),rownames(MIA_vs_AIS))
##有一些重复的差异基因你
DEG = b[rownames(b) %in% DEGgene,]
str(DEG)
dat = as.data.frame(t(DEG))
pca.res <- prcomp(dat, scale. = T, center = T)
pca.res
tmp <- as.data.frame(pca.res$x)
head(tmp)

library(scatterplot3d)
rownames(tmp)
col = ifelse(str_detect(rownames(tmp),"AIS"),'purple',
             ifelse(str_detect(rownames(tmp),"IAC"),'orange','green'))
scatterplot3d(tmp[,1:3],color=col,
              pch = 16,angle=30,
              box=T,type="p",
              lty.hide=2,lty.grid = 2)
legend("topleft",c('AIS','IAC','MIA'),
       fill=c('purple','orange','green'),box.col=NA,cex=0.7)

###文章中PCA显示MIA和IAC在转录组水平上非常接近
##我这里显示的是IAC和AIS更为接近，倒是和前面的韦恩图对应上了
##但是这更能说明IAC是最终状态啊，期间的MIA是中间状态，感觉也是可以解释的




###两组差异分析共显示dotplot图####
#假设我们有某个癌症的组织样本，cancer vs normal得到癌变后差异表达蛋白，而同样的。我们也去取血液的样本，检测血液中蛋白的变化，找到差异蛋白，这个时候就可以将这两个数据合在一起比较，比如某些基因在癌组织和血液中变化一致，可以认为他们是正常相关的，通过这样的思路我们可以确定某些标志物，也就是血液能够代替活检的相关蛋白筛选。
#IAC_vs_MIA,MIA_vs_AIS
table(Idents(sce))
as.data.frame(table(Idents(sce)))

IACMIA = all_markers_sig[[1]]
MIAAIS = all_markers_sig[[3]]

ids=intersect(rownames(IACMIA),
              rownames(MIAAIS))
library(dplyr)
library(magrittr)
###数量太多了，随机选取100个，不具有实际意义，可以按照绝对值来挑选
#dat = IACMIA[ids,]
#top_100 <- dat %>% arrange(desc(abs(avg_log2FC))) %>% head(100)
id = sample(ids,200)

###这里又出现了另一个问题，要按照每个细胞亚群的差异基因来分组group,6个主要亚群
sce.markers = read.csv('../4-plot/sce.markers.csv',row.names = 1)
markergene = sce.markers[sce.markers$gene %in% id,]
##有的基因是不同的细胞亚群中是重复的(注意：可以自己根据实际意义来选择去掉哪一项，或者删掉重复基因)
markergene = markergene[!duplicated(markergene$gene),]
df= data.frame(
  IACMIA = IACMIA[markergene$gene,'avg_log2FC'],
  MIAAIS = MIAAIS[markergene$gene,'avg_log2FC'],
  gene = markergene$gene,
  celltype = markergene$cluster
)


library(ggpubr)
ggscatter(df, x = "MIAAIS", y = "IACMIA",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson",  label.sep = "\n")
)

library(ggplot2)
library(ggrepel)
library(ggthemes)
A=df
ggplot(A, aes(x=MIAAIS, y=IACMIA,color = celltype)) +
  geom_hline(yintercept= c(0, 0), color = "black",  size=1) +#添加横线
  geom_vline(xintercept=c(0, 0), color = "black", size=1)+
  geom_point(size = 3,shape=21)

colnames(A)
ggplot(A, aes(x=MIAAIS, y=IACMIA,color = celltype)) +
  geom_hline(yintercept= c(0, 0), color = "black",  size=0.5) +#添加横线
  geom_vline(xintercept=c(0, 0), color = "black", size=0.5)+
  geom_point(size = 3)+
  xlim(-2,3)+
  ylim(-3, 2)+
  labs(x = "Log2FC MIA enriched(MIA vs AIS))",
       y = "Log2FC IAC enriched(IAC vs MIA)", title = "") + 
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5), 
        axis.text = element_text(size = 14, color = 'black'), 
        axis.title = element_text(size = 14, color = 'black'))+
  theme(legend.position = "right")+
  theme_few()+
  geom_text_repel(data=A, aes(label=gene), color="black", size=2, fontface="italic", 
                  point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)


 
###三组差异分析火山图####
library(tidyverse)

all_markers_sig2 = lapply(markers_list, function(x){
  markers_sig <- subset(x,abs(avg_log2FC)>0.5)
})

IAC_vs_MIA2 = all_markers_sig2[[1]]
IAC_vs_AIS2 = all_markers_sig2[[2]]
MIA_vs_AIS2 = all_markers_sig2[[3]]
IAC_vs_MIA2$group = 'IAC_vs_MIA2'
IAC_vs_MIA2$gene = rownames(IAC_vs_MIA2)
MIA_vs_AIS2$group = 'MIA_vs_AIS2'
MIA_vs_AIS2$gene = rownames(MIA_vs_AIS2)
IAC_vs_AIS2$group = 'IAC_vs_AIS2'
IAC_vs_AIS2$gene = rownames(IAC_vs_AIS2)
dat = rbind(MIA_vs_AIS2,IAC_vs_AIS2,IAC_vs_MIA2)

#添加显著性标签：
colnames(dat)
dat$label <- ifelse(dat$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
head(dat)
table(dat$label)
#依次获取最显著的基因
top_15_rows <- dat %>%
  group_by(group) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n = 15)
table(top_15_rows$group)

#新增一列，将Top差异基因标记为2，其他的标记为1
dat$size <- case_when(!(dat$gene %in% top_15_rows$gene)~ 1,
                      dat$gene %in% top_15_rows$gene ~ 2)
table(dat$size)
#提取非Top10的基因表格；
dt <- filter(dat,size==1)
#绘制散点火山图
dt <- filter(dat,size==1)
head(dt)
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = group, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)
p


#叠加每个Cluster Top10基因散点
str(dt)
dt = as.data.frame(dt)
p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = group, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top_15_rows,
              aes(x = group, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)
p

#根据图p中log2FC区间确定背景柱长度
dfbar<-data.frame(x=c(1,2,3),
                  y=c(10,10,11))
dfbar1<-data.frame(x=c(1,2,3),
                   y=c(-7.5,-7,-5))
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

#把散点火山图叠加到背景柱上
p2 <- ggplot() +
  geom_col(data = dfbar, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = dfbar1, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_jitter(data = dt, aes(x = group, y = avg_log2FC, color = label), size = 0.85, width = 0.4) +
  geom_jitter(data = top_15_rows, aes(x = group, y = avg_log2FC, color = label), size = 1, width = 0.4) +
  scale_x_discrete() 
p2

#添加X轴的stage色块标签
dfcol<-data.frame(x=c(1:3),
                  y=0,
                  label=c('MIA_vs_AIS','IAC_vs_AIS','IAC_vs_MIA'))
mycol <- c("#00A0877F","#3C54887F","#F39B7F7F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 1.7,
                     show.legend = F)
p3


#给每个stage差异表达前Top15基因加上标签
p4 <- p3+
  geom_text_repel(
    data=top_15_rows,
    aes(x = group, y = avg_log2FC,label = gene),
    force = 1.2,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last")
  )
p4

#散点颜色调整
p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p5


#修改X/Y轴标题和添加cluster数字：
p6 <- p5+
  labs(x="Cluster",y="average logFC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6

#自定义主题美化：
p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7


setwd('../')
