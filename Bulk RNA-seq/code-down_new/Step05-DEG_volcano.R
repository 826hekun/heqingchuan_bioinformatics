rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(tidyverse)
 
# 读差异分析结果
lname <- load(file = "data//Step03-edgeR_nrDEG.Rdata")

# 根据需要修改DEG的值
data <- DEG_edgeR_symbol
colnames(data)


# 绘制火山图
colnames(data)
p <- ggplot(data=data, aes(x=logFC, y=-log10(PValue),color=regulated)) + 
  geom_point(alpha=0.5, size=1.2) + 
  theme_set(theme_set(theme_bw(base_size=20))) + theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  xlab("log2FC") + ylab("-log10(Pvalue)") +
  scale_colour_manual(values = c(down='blue',normal='grey',up='red')) +
  geom_vline(xintercept=c(-(log2(1.5)),log2(1.5)),lty=2,col="black",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.6)
p


# 添加top基因
# 通过FC选取TOP10
label <- data[order(abs(data$logFC),decreasing = T)[1:10],]
# 通过pvalue选取TOP10
#label <- data[order(abs(data$PValue),decreasing = F)[1:10],]
label <- na.omit(label)
label

p1 <- p + geom_point(size = 3, shape = 1, data = label) +
  ggrepel::geom_text_repel( aes(label = SYMBOL), data = label, color="black" )

p1


# 保存结果
png(file = "result/5.Volcano_Plot.png",width = 900, height = 800, res=150)
plot(p1)
dev.off()




