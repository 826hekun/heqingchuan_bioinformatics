options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

cran_packages <- c('tidyr',
                   'tibble',
                   'dplyr',
                   'stringr',
                   'ggplot2',
                   'ggpubr',
                   'factoextra',
                   'FactoMineR',
                   'pheatmap',
                   "survival",
                   "survminer",
                   "patchwork",
                   "statsExpressions",
                   "ggstatsplot",
                   "ggplotify",
                   "cowplot",
                   "glmnet",
                   "ROCR",
                   "caret",
                   "randomForest",
                   "survminer",
                   "Hmisc",
                   "e1071",
                   "deconstructSigs",
                   "AnnoProbe",
                   "timeROC",
                   "circlize",
                   "VennDiagram",
                   "tinyarray"
) 
Biocductor_packages <- c("limma",
                         "clusterProfiler",
                         "org.Hs.eg.db",
                         "SummarizedExperiment",
                         "DESeq2",
                         "edgeR",
                         "ggpubr",
                         "rtracklayer",
                         "genefilter",
                         "maftools",
                         "ComplexHeatmap",
                         "TCGAbiolinks"
)


for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


# use BiocManager to install
for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

#前面的任何信息都先不要管。主要看这里
for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}

#没有任何提示就是成功了，如果有warningxx包不存在，用library检查一下。

#报错就回去重新安装。如果你没有安装xx包，却提示你xx包不存在，这也正常，是因为依赖关系，缺啥补啥。

