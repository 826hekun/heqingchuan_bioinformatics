rm(list=ls())
###根据TCGA-LUAD验证细胞亚群的预后生存
dir.create("11-TCGA_LUAD")
setwd('11-TCGA_LUAD/')
source('../scRNA_scripts/mycolors.R')
library(tidyverse)
library(data.table) 
#如果需要下载TCGA其他癌种，请更换proj
proj = "TCGA-LUAD"
dir.create("input")

##由于网络限制，这里我下载好了数据存放在input文件夹
##如果需要自己下载，把F改为T
if(F){
  download.file(url = paste0("https://gdc.xenahubs.net/download/",proj, ".htseq_counts.tsv.gz"),destfile = paste0("input/",proj,".htseq_counts.tsv.gz"))  ##表达数据
  download.file(url = paste0("https://gdc.xenahubs.net/download/",proj, ".GDC_phenotype.tsv.gz"),destfile = paste0("input/",proj,".GDC_phenotype.tsv.gz")) ##临床数据
  download.file(url = paste0("https://gdc.xenahubs.net/download/",proj, ".survival.tsv"),destfile = paste0("input/",proj,".survival.tsv")) ##生存数据
}

clinical = read.delim(paste0("input/",proj,".GDC_phenotype.tsv.gz"),fill = T,header = T,sep = "\t")
surv = read.delim(paste0("input/",proj,".survival.tsv"),header = T) 
head(surv) #生存数据os和os.time

### 1.处理表达矩阵和分组信息
#### 1.1 表达矩阵
library(data.table) 
dat <- data.table::fread('input/TCGA-LUAD.htseq_counts.tsv.gz',
                         data.table = F)  
head(dat[,1:4])
tail(dat[,1:4]) 
dat = dat[1:(nrow(dat)-5),]
rownames(dat) = dat$Ensembl_ID
a = dat
a = a[,-1]
##逆转 log
a = as.matrix(2^a - 1)
# 用apply转换为整数矩阵
head(a[,1:4])
tail(a[,1:4]) 
colSums(a)/1e6
exp = apply(a, 2, as.integer)
rownames(exp) = rownames(dat)
exp= log(edgeR::cpm(exp)+1)
library(stringr)
head(rownames(exp))
library(AnnoProbe)
library(tinyarray)
rownames(exp) = substr(rownames(exp), 1, 15) 
re = annoGene(rownames(exp),ID_type = "ENSEMBL");head(re)
exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")
head(exp[,1:4])
tail(exp[,1:4]) 
proj='tcga-luad'
save(exp,file = paste0(proj,".htseq_counts.rdata") ) 



rm(list=ls())
proj='tcga-luad'
load(file = paste0(proj,".htseq_counts.rdata") )
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')  
Group = factor(Group,levels = c("normal","tumor"))
print(table(Group))
# 生存分析只需要tumor样品即可
exprSet = exp[,Group=='tumor']

clinical = read.delim('input/TCGA-LUAD.GDC_phenotype.tsv.gz',
                      fill = T,header = T,sep = "\t")
surv = read.delim('input/TCGA-LUAD.survival.tsv',header = T) 
library(tidyverse)
meta = left_join(surv,clinical,by = c("sample"= "submitter_id.samples"))
head(meta[,1:4])
tail(meta[,1:4]) 
print(dim(meta))

#去掉生存信息不全或者生存时间小于30天的样本，样本纳排标准不唯一，且差别很大.
k1 = meta$OS.time >= 30
k2 = !(is.na(meta$OS.time)|is.na(meta$OS))
meta = meta[k1&k2,]
meta = meta[,c(
  'sample',
  'OS',
  'OS.time'
)]
colnames(meta)=c('ID','event','time')
meta$time = meta$time/30
rownames(meta) <- meta$ID
s = intersect(rownames(meta),colnames(exprSet))
exprSet = exprSet[,s]
meta = meta[s,]
identical(rownames(meta),colnames(exprSet))
save(exprSet,meta,file = paste0(proj,".for_survival.rdata") ) 


##基于单细胞转录组的生存分析
#第1步：根据单细胞亚群基因集在肿瘤病人表达量矩阵里面进行gsva打分
rm(list=ls())
library(survival)
library(survminer) 
library(ggstatsplot) 
library(gplots)
library(ggplot2) 
library(pheatmap)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(GSVA) 
library(GSEABase)

# 1. 载入表达量矩阵和临床信息 ----
proj='tcga-luad'
load(file = paste0(proj,".for_survival.rdata") ) 
exprSet[1:4,1:4]
phe=meta
head(phe)
mySurv <- with(phe, Surv(time, event))
survival_dat=phe

# 2. creat geneset----
markers = read.csv('../7-epi/epimarkers.csv',row.names = 1)
head(markers$cluster)
markers$cluster <- as.factor(markers$cluster)

deg_list <- split(markers$gene, markers$cluster)
deg_list

gs = lapply(deg_list, toupper) 
geneset <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, gs, names(gs)))
geneset

# 3. run gsva----
X=as.matrix(exprSet)
es.max <- gsva(X, geneset, 
               mx.diff=FALSE, verbose=FALSE, 
               parallel.sz=4)
es.max[1:4, 1:4] 
pheatmap(es.max) 


#第2步：根据gsva打分值高低分组进行生存分析
#一般来说， 就是根据gsva打分的中位值高低分组进行生存分析，代码如下所示：

# 4. 根据gsva结果高低分组后批量生存分析 ----
es.max[1:4, 1:4]
splots <- list()
g = 1
for (i in  names(deg_list) ) {
  # i =  names(deg_list) [1]
  subset = paste0('cluster_',i)
  print(subset)
  v = as.numeric(es.max[i,])   #每一个亚群表达量。
  sub_group <- ifelse( v < 0,"low","high")   #如果表达量小于0的话，就定义为low。gsva处理过表达量。0.几左右
  table(sub_group) 
  phe$sub_group=sub_group
  # Fit survival curves
  require("survival")
  fit <- survfit(Surv(time, event) ~ sub_group, data = phe)
  library("survminer")
  survp <- ggsurvplot(fit, data = phe,
                      surv.median.line = "hv", # Add medians survival
                      pval = TRUE,             # Add p-value and tervals 
                      conf.int = TRUE,        # Add the 95% confidence band
                      risk.table = TRUE,      # Add risk table
                      tables.height = 0.2,
                      tables.theme = theme_cleantable(),
                      palette = "jco",
                      ggtheme = theme_bw(),
                      title = subset)
  print(survp)
  splots[[g]] <-  survp
  g = g + 1
}

length(splots)
x1 = ceiling(sqrt(length(splots)))
y1 = x1

all_plot <- arrange_ggsurvplots(splots,
                                print = F,
                                ncol = x1, nrow = y1,
                                risk.table.height = 0.3,
                                surv.plot.height = 0.7)
all_plot 
x2=5*x1
y2=5*y1
prefix=''
pro=''
ggsave(all_plot, #path = prefix,
       filename = paste0(pro, 'all_survival_plot.pdf'),
       width = x2,height = y2)

#就可以看到如下所示的每个单细胞亚群的生存分析结果，很明显跟文章类似的，也是增殖亚群（UBE2C+） 它这个亚群的基因可以在tcga数据库的luad数据集里面的有统计学显著的生存分析意义。


#第3步：根据gsva打分值进行取巧分组进行生存分析
#如果是上面的根据gsva打分的中位值高低分组进行生存分析都没有生存分析统计学显著意义，但是又想看看每个亚群的具体的到底是保护因子还是风险因子，也可以使用surv_cutpoint函数哦：
## 感觉这种刻意找显著差异的做法maybe会受到质疑，但是存在即合理
head(phe)
csplots <- list()
cg = 1
for (i in  names(deg_list) ) {
  # i =  names(deg_list) [1]
  subset = paste0('cluster_',i)
  print(subset)
  v = as.numeric(es.max[i,])   #每一个亚群表达量。
  phe$v <- v
  head(phe)
  sur.cut <- surv_cutpoint(phe,
                           time= 'time',
                           event = 'event' ,
                           variables = 'v' )
  sur.cat <- surv_categorize(sur.cut)
  head(sur.cat)
  sfit <- survfit(Surv(time, event) ~ v, data = sur.cat)
  p_surv_cut <- ggsurvplot(sfit, data = phe,
                           surv.median.line = "hv", # Add medians survival
                           pval = TRUE,             # Add p-value and tervals 
                           conf.int = TRUE,        # Add the 95% confidence band
                           risk.table = TRUE,      # Add risk table
                           tables.height = 0.2,
                           tables.theme = theme_cleantable(),
                           palette = "jco",
                           ggtheme = theme_bw(),
                           title = subset)
  print(p_surv_cut)
  csplots[[cg]] <-  p_surv_cut
  cg = cg + 1
}

length(csplots)
x1 = ceiling(sqrt(length(csplots)))
y1 = x1

all_plot <- arrange_ggsurvplots(csplots,
                                print = F,
                                ncol = x1, nrow = y1,
                                risk.table.height = 0.3,
                                surv.plot.height = 0.7)
all_plot 
x2=5*x1
y2=5*y1
ggsave(all_plot, #path = prefix,
       filename = paste0(pro, 'all_cut_point_survival_plot.pdf'),width = x2,height = y2)

#因为它并不是中位值高低分组，所以两个分组的病人数量是不平衡的




setwd('../')



