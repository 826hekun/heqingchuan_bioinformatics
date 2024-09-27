# 魔幻操作，一键清空
rm(list = ls()) 
options(stringsAsFactors = F)

library(FactoMineR)
library(factoextra)
library(corrplot)
library(pheatmap)
library(tidyverse)
 
# 加载数据并检查
lname <- load(file = 'data/Step01-airwayData.Rdata')
lname



## 1.样本之间的相关性-层次聚类树----
dat <- log10(express_cpm+1)
dat[1:4,1:4]
dim(dat)
sampleTree <- hclust(dist(t(dat)), method = "average")
plot(sampleTree)

# 提取样本聚类信息
temp <- as.data.frame(cutree(sampleTree,k = 2)) %>% 
  rownames_to_column(var="sample")
temp1 <- merge(temp,group,by.x = "sample",by.y="run_accession")
table(temp1$`cutree(sampleTree, k = 2)`,temp1$sample_title)

# 保存结果
pdf(file = "result/2.sample_Treeplot.pdf",width = 7,height = 6)
plot(sampleTree)
dev.off()



## 2.样本之间的相关性-PCA----
# 第一步，数据预处理
dat <- log10(express_cpm+1)
dat[1:4,1:4]

dat <- as.data.frame(t(dat))
dat_pca <- PCA(dat, graph = FALSE)

group_list <- group[match(group$run_accession,rownames(dat)), 2]
group_list

# geom.ind: point显示点，text显示文字
# palette: 用不同颜色表示分组
# addEllipses: 是否圈起来
mythe <- theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

p <- fviz_pca_ind(dat_pca,
                  geom.ind = "text",  #text
                  col.ind = group_list, 
                  palette = c("#00AFBB", "#E7B800"), 
                  addEllipses = T,  
                  legend.title = "Groups") + mythe
p

# 保存结果
pdf(file = "result/2.sample_PCA.pdf",width = 6.5,height = 6)
plot(p)
dev.off()




## 3.样本之间的相关性-cor----
# 选择差异变化大的基因算样本相关性
exprSet <- express_cpm
exprSet = exprSet[names(sort(apply(exprSet, 1, mad),decreasing = T)[1:800]),]
dim(exprSet)

# 计算相关性
M <- cor(exprSet,method = "spearman")
M

# 构造注释条
anno <- data.frame(group=group$sample_title,row.names = group$run_accession )

# 保存结果
pheatmap(M,display_numbers = T,
         annotation_col = anno,
         fontsize = 10,cellheight = 30,
         cellwidth = 30,cluster_rows = T,
         cluster_cols = T,
         filename = "result/2.sample_Cor.pdf",width = 7.5,height = 7)




