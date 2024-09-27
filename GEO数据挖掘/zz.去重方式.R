rm(list = ls())
load("step2output.Rdata")
#保留最大值
exp2 = exp[ids$probe_id,]
identical(ids$probe_id,rownames(exp2))
library(dplyr)
ids = ids %>% 
  mutate(exprowsum = rowSums(exp2)) %>% 
  arrange(desc(exprowsum)) %>% 
  select(-3) %>% 
  distinct(symbol,.keep_all = T)
nrow(ids)
# 拿这个ids去inner_join

#求平均值
rm(list = ls())
load("step2output.Rdata")
exp3 = exp[ids$probe_id,]
rownames(exp3) = ids$symbol
exp3[1:4,1:4]
exp4 = limma::avereps(exp3)

# 此时拿到的exp4已经是一个基因为行名的表达矩阵，直接差异分析，不再需要inner_join 
