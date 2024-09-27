####构建左下角坐标轴
left_axes = function(scRNA){
# extact PC ranges
pc12 <- Embeddings(object = scRNA,reduction = 'umap') %>%
  data.frame()

# check
head(pc12,3)
# get botomn-left coord
lower <- floor(min(min(pc12$umap_1),min(pc12$umap_2))) - 2

# get relative line length
linelen <- abs(0.3*lower) + lower

# mid point
mid <- abs(0.3*lower)/2 + lower

# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('umap_2','umap_1'),each = 2))

# axies label
label <- data.frame(lab = c('umap_2','umap_1'),angle = c(90,0),
                    x = c(lower - 3,mid),y = c(mid,lower - 2.5))

# 组合坐标轴和标签数据
result <- list(axes = axes, label = label)

# 返回组合后的数据框
return(result)
}

