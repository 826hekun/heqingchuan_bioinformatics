rm(list = ls())
library(stringr)
 
## ====================1.读取数据
# 读取raw count表达矩阵
rawcount <- read.table("data/raw_counts.txt",row.names = 1, 
                       sep = "\t", header = T)
colnames(rawcount)

# 查看表达谱
rawcount[1:6, 1:6]

# 去除前的基因表达矩阵情况
dim(rawcount)

# 获取分组信息
# 差异分析方案为：Dex vs untreated
group <- read.table("data/group.txt", header = T,sep = "\t", quote = "\"")
group


## =================== 2.表达矩阵预处理
# 过滤低表达基因
keep <- rowSums(rawcount>0) >= floor(0.75*ncol(rawcount))
table(keep)

filter_count <- rawcount[keep,]
filter_count[1:4,1:4]
dim(filter_count)

# 加载edgeR包计算counts per millio(cpm) 表达矩阵
library(edgeR)
express_cpm <- cpm(filter_count)
express_cpm[1:6,1:6]

# 保存表达矩阵和分组结果
save(filter_count, express_cpm, group, file = "data/Step01-airwayData.Rdata")


