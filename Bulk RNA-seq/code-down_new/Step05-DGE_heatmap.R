rm(list = ls())
options(stringsAsFactors = F)
# 加载包
library(pheatmap) #complexheatmap
library(tidyverse)
 
# 加载原始表达矩阵
lname <- load(file = "data/Step01-airwayData.Rdata")
lname

express_cpm1 <- rownames_to_column(as.data.frame(express_cpm) ,var = "ID")

# 读取差异分析结果
lname <- load(file = "data/Step03-edgeR_nrDEG.Rdata")
lname

# 提取所有差异表达的基因名
edgeR_sigGene <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated!="normal",]
head(edgeR_sigGene)

data <- merge(edgeR_sigGene,express_cpm1,by.x = "ENSEMBL",by.y = "ID")
data <- na.omit(data)
data <- data[!duplicated(data$SYMBOL),]

# 绘制热图
dat <- select(data,starts_with("SRR"))
rownames(dat) <- data$SYMBOL
dat[1:4,1:4]
anno <- data.frame(group=group$sample_title,row.names = group$run_accession)

pheatmap(dat,scale = "row",show_colnames =T,
         show_rownames = F, cluster_cols = T,
         annotation_col=anno,
         main = "edgeR's DEG")


# 显示指定symbol,这里随便展示10个基因symbol
labels <- rep(x = "",times=nrow(dat))
labels[1:10] <- rownames(dat)[1:10]
pheatmap(dat,scale = "row",show_colnames =T,
         show_rownames = T, 
         cluster_cols = T,
         annotation_col=anno,
         labels_row = labels,
         fontsize_row = 8,
         main = "edgeR's DEG")



# 按照指定顺序绘制热图
dex_exp <- express_cpm[,match(rownames(anno)[which(anno$group=="Dex")],
                              colnames(express_cpm))]

untreated_exp <- express_cpm[,match(rownames(anno)[which(anno$group=="untreated")],
                              colnames(express_cpm))]

data_new <- cbind(dex_exp, untreated_exp)
dat1 <- data_new[match(edgeR_sigGene$ENSEMBL,rownames(data_new)),]

pheatmap(dat1, scale = "row",show_colnames =T,show_rownames = F, 
              cluster_cols = F, 
              annotation_col=anno,
              main = "edgeR's DEG")

