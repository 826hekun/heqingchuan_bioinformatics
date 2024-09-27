rm(list = ls())
options(stringsAsFactors = F)

# 加载包，设置绘图参数
library(ggplot2)
library(ggsci)
library(tidyverse)

mythe <- theme_bw() + theme(panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank())


# 加载原始表达的数据
lname <- load(file = "data/Step01-airwayData.Rdata")
lname
exprSet <- log10(as.matrix(express_cpm)+1)
exprSet[1:6,1:6]



## 1.样本表达总体分布-箱式图
# 构造绘图数据
data <- exprSet %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = "sample",values_to = "expression")

head(data)

p <- ggplot(data = data, aes(x=sample,y=expression,fill=sample))
p1 <- p + geom_boxplot() + 
  mythe + theme(axis.text.x = element_text(angle = 90)) + 
  xlab(NULL) + ylab("log10(CPM+1)") + scale_fill_lancet()

p1

# 保存图片
png(file = "result/1.sample_boxplot.png",width = 800, height = 900,res=150)
print(p1)
dev.off()



## 2.样本表达总体分布-小提琴图
p2 <- p + geom_violin() +  mythe +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90)) + 
  xlab(NULL) + ylab("log10(CPM+1)")+scale_fill_lancet()
p2

# 保存图片
png(file = "result/1.sample_violin.png",width = 800, height = 900,res=150)
print(p2)
dev.off()



## 3.样本表达总体分布-概率密度分布图
m <- ggplot(data=data, aes(x=expression)) 
p3 <- m +  geom_density(aes(fill=sample, colour=sample),alpha = 0.1) + 
  xlab("log10(CPM+1)") + mythe + scale_fill_npg()
p3

# 保存图片
png(file = "result/1.sample_density.png",width = 800, height = 700, res=150)
print(p3)
dev.off()

