#设置镜像
options("repos"=c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#安装R包
if(!require(ggplot2))install.packages('ggplot2',update = F,ask = F)
if(!require(ggpubr))install.packages('ggpubr',update = F,ask = F)
if(!require(eoffice))install.packages("eoffice",update = F,ask = F)
if(!require(patchwork))install.packages("patchwork",update = F,ask = F)
#加载以检查是否安装成功
library(ggplot2)
library(ggpubr)
library(eoffice)
library(patchwork)
