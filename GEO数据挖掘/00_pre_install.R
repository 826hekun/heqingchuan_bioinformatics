options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

cran_packages <- c('tidyr',
                   'tibble',
                   'dplyr',
                   'stringr',
                   'ggplot2',
                   'ggpubr',
                   'factoextra',
                   'FactoMineR',
                   'devtools',
                   'cowplot',
                   'patchwork',
                   'basetheme',
                   'paletteer',
                   'AnnoProbe',
                   'ggthemes',
                   'VennDiagram') 
Biocductor_packages <- c('GEOquery',
                         'hgu133plus2.db',
                         'ggnewscale',
                         "limma",
                         "impute",
                         "GSEABase",
                         "GSVA",
                         "clusterProfiler",
                         "org.Hs.eg.db",
                         "preprocessCore",
                         "enrichplot")

for (pkg in cran_packages){
  if (! require(pkg,character.only=T,quietly = T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T,quietly = T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

#前面的所有提示和报错都先不要管。主要看这里
for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}
#没有任何提示就是成功了，如果有warning xx包不存在，用library检查一下。
#library报错，就单独安装。
if(!require(tinyarray))install.packages("https://cran.r-project.org/src/contrib/Archive/tinyarray/tinyarray_2.3.3.tar.gz",repos = NULL)
library(tinyarray)
