jimmy <- function(a,b,m = 2){
  (a+b)^m
}
jimmy(a = 1,b = 2)
jimmy(1,2)
jimmy(3,6)
jimmy(3,6,-2)

#复习：绘图函数plot()
par(mfrow = c(2,2)) #把画板分成四块，两行两列
#如果报错，把右下角画板拉大一点即可
x = c(2,5,6,2,9);plot(x)
x = seq(2,80,4);plot(x)
x = rnorm(10);plot(x)
x = iris$Sepal.Length;plot(x)

#思考：plot画iris的前四列？
plot(iris[,1],col = iris[,5])
plot(iris[,2],col = iris[,5])
plot(iris[,3],col = iris[,5])
plot(iris[,4],col = iris[,5])

#当一个代码需要复制粘贴三次，就应该写成函数或使用循环

jimmy <- function(i){
  plot(iris[,i],col=iris[,5])
}

jimmy(1)
jimmy(2)
jimmy(3)
jimmy(4)

# R包安装

options("repos"=c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")

install.packages("tidyr")
install.packages('BiocManager')
BiocManager::install("ggplot2")
install.packages('devtools')
devtools::install_github("jmzeng1314/idmap1") #括号里写作者用户名加包名

# 清华镜像
# http://mirrors.tuna.tsinghua.edu.cn/CRAN/
# http://mirrors.tuna.tsinghua.edu.cn/bioconductor/
  
# 中科大镜像
# http://mirrors.ustc.edu.cn/CRAN/
# http://mirrors.ustc.edu.cn/bioc/
  
library(tidyr)
require(tidyr)

# 分情况讨论

if(!require(stringr))install.packages("stringr")

# 获取帮助
?seq
library(stringr)
browseVignettes("stringr")
ls("package:stringr")

