rm(list = ls())
if(!require(stringr))install.packages('stringr')
library(stringr)

x <- "The birch canoe slid on the smooth planks."
x
### 1.检测字符串长度
str_length(x)
length(x)
### 2.字符串拆分
str_split(x," ")
x2 = str_split(x," ")[[1]];x2

y = c("jimmy 150","nicker 140","tony 152")
str_split(y," ")
str_split(y," ",simplify = T)

### 3.按位置提取字符串
str_sub(x,5,9)

### 4.字符检测
str_detect(x2,"h")
str_starts(x2,"T")
str_ends(x2,"e")
### 5.字符串替换
x2
str_replace(x2,"o","A")
str_replace_all(x2,"o","A")

### 6.字符删除
x
str_remove(x," ")
str_remove_all(x," ")
