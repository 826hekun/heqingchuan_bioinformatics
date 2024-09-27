rm(list = ls())
options(stringsAsFactors = F)

# 加载包
library(edgeR)
library(ggplot2)
 
# 读取基因表达矩阵信息并查看分组信息和表达矩阵数据
lname <- load(file = "data/Step01-airwayData.Rdata")
lname

# 表达谱
filter_count[1:4,1:4]

# 分组信息
group_list <- group[match(colnames(filter_count),group$run_accession),2]
group_list

# treat vs control
comp <- unlist(strsplit("Dex_vs_untreated",split = "_vs_"))
group_list <- factor(group_list,levels = comp)
group_list
table(group_list)


# 构建线性模型。0代表x线性模型的截距为0
design <- model.matrix(~0+group_list)
rownames(design) <- colnames(filter_count)
colnames(design) <- levels(factor(group_list))
design

# 构建edgeR的DGEList对象
DEG <- DGEList(counts=filter_count, 
               group=factor(group_list))

# 归一化基因表达分布
DEG <- calcNormFactors(DEG)

# 计算线性模型的参数
DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

# 拟合线性模型
fit <- glmFit(DEG, design)

# 进行差异分析
lrt <- glmLRT(fit, contrast=c(1,-1)) 

# 提取过滤差异分析结果
DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG),adjust.method = "BH"))
head(DEG_edgeR)

# 筛选上下调，设定阈值
fc_cutoff <- 1.5 # 1.2  1.5  2
pvalue <- 0.05  # 0.05, 0.01, 0.001
#fdr <- 0.05  # 0.05, 0.01, 0.001

DEG_edgeR$regulated <- "normal"

loc_up <- intersect(which( DEG_edgeR$logFC > log2(fc_cutoff) ),
                    which( DEG_edgeR$PValue < pvalue) )

loc_down <- intersect(which(DEG_edgeR$logFC < (-log2(fc_cutoff))),
                      which(DEG_edgeR$PValue<pvalue))

DEG_edgeR$regulated[loc_up] <- "up"
DEG_edgeR$regulated[loc_down] <- "down"

table(DEG_edgeR$regulated)


## 添加一列gene symbol
# 方法1：使用包
library(org.Hs.eg.db) #不同物种 ，包不一样
keytypes(org.Hs.eg.db)

library(clusterProfiler)
id2symbol <- bitr(rownames(DEG_edgeR), 
                  fromType = "ENSEMBL", 
                  toType = "SYMBOL", 
                  OrgDb = org.Hs.eg.db)
head(id2symbol)

DEG_edgeR <- cbind(GeneID=rownames(DEG_edgeR),DEG_edgeR)
DEG_edgeR_symbol <- merge(id2symbol,DEG_edgeR,
                          by.x="ENSEMBL",by.y="GeneID",all.y=T)
head(DEG_edgeR_symbol)


# 方法2：gtf文件中得到的id与name关系
# Assembly: GRCh37(hg19) Release: ？
# 使用上课测试得到的count做



# 选择显著差异表达的结果
library(tidyverse)
DEG_edgeR_symbol_Sig <- filter(DEG_edgeR_symbol,regulated!="normal")

# 保存
write.csv(DEG_edgeR_symbol,"result/4.DEG_edgeR_all.csv", row.names = F)
write.csv(DEG_edgeR_symbol_Sig,"result/4.DEG_edgeR_Sig.csv", row.names = F)
save(DEG_edgeR_symbol,file = "data/Step03-edgeR_nrDEG.Rdata")



##====== 检查是否上下调设置错了
# 挑选一个差异表达基因
head(DEG_edgeR_symbol_Sig)

exp <- c(t(express_cpm[match("ENSG00000001626",rownames(express_cpm)),]))
test <- data.frame(value=exp, group=group_list)

ggplot(data=test,aes(x=group,y=value,fill=group)) + geom_boxplot()


