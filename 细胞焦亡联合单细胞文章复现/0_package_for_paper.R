options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor") 

cran_packages <- c("Seurat","harmony","tidyverse","ggplot2",
                   "paletteer","ggpubr","rms","oncoPredict",
                   "gtools","reshape2",
                   'factoextra','FactoMineR',
                   'pheatmap',"survival","survminer",
                   "patchwork","statsExpressions",
                   "ggstatsplot", "ggplotify",
                   "cowplot","glmnet","ROCR",
                   "caret","randomForest","Hmisc",
                   "e1071","deconstructSigs",
                   "AnnoProbe","timeROC",
                   "circlize","VennDiagram",
                   "tinyarray") 
Biocductor_packages <- c("celldex","SingleR",
                         "BiocParallel","GSEABase",
                         "AUCell","maftools","limma",
                         "clusterProfiler","org.Hs.eg.db",
                         "SummarizedExperiment","DESeq2",
                         "edgeR","ggpubr",
                         "rtracklayer","genefilter",
                         "ComplexHeatmap","TCGAbiolinks")


for (pkg in cran_packages){
  if (! require(pkg,character.only=T,quietly = T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


# use BiocManager to install
for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T,quietly = T) ) {
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

if(!require(CIBERSORT))devtools:: install_github ("Moonerss/CIBERSORT")
if(!require(TCGAmutations))devtools::install_github("PoisonAlien/TCGAmutations")
if(!require(tinyarray))install.packages("https://cran.r-project.org/src/contrib/Archive/tinyarray/tinyarray_2.3.3.tar.gz",repos = NULL)
