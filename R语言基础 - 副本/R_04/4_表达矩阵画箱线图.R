# 表达矩阵
set.seed(10086)
exp = matrix(rnorm(18),ncol = 6)
exp = round(exp,2)
rownames(exp) = paste0("gene",1:3)
colnames(exp) = paste0("test",1:6)
exp[,1:3] = exp[,1:3]+1
exp

library(tidyr)
library(tibble)
library(dplyr)
dat = t(exp) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group = rep(c("control","treat"),each = 3))

pdat = dat%>% 
  pivot_longer(cols = starts_with("gene"),
               names_to = "gene",
               values_to = "count")

library(ggplot2)
p = ggplot(pdat,aes(gene,count))+
  geom_boxplot(aes(fill = group))+
  theme_bw()
p
p + facet_wrap(~gene,scales = "free")



