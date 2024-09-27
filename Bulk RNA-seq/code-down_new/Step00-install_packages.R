rm(list = ls())

## Installing R packages
bioPackages <-c( "corrplot", "ggrepel", "stringr", "FactoMineR",
  "factoextra", "limma", "pheatmap", "edgeR", "DESeq2", "clusterProfiler",
  "org.Hs.eg.db", "GSEABase", "tidyverse", "GSVA" )

## If you are in China, run the command below
#options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/")) 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.bfsu.edu.cn/CRAN/")) 

# cannot open URL 'https://mirrors.ustc.edu.cn/... 时修改为这个
options(download.file.method = 'libcurl') 
options(url.method='libcurl')

# 检查是否设定完毕
options()$repos 
options()$BioC_mirror

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 安装devtools管理github上的软件包
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")


## Installing missing packages
lapply( bioPackages, function( bioPackage ){
  if(!bioPackage %in% rownames(installed.packages())){
    CRANpackages <- available.packages()
    if(bioPackage %in% rownames(CRANpackages)){
      install.packages( bioPackage)
    }else{
      BiocManager::install(bioPackage,suppressUpdates=F,ask=F)
    }
  }
})


## 验证R包是否安装成功
library(limma)
library(edgeR)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(clusterProfiler)
library(org.Hs.eg.db)

