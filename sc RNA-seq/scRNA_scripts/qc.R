basic_qc <- function(input_sce){
  #计算线粒体基因比例
  mito_genes=rownames(input_sce)[grep("^MT-", rownames(input_sce),ignore.case = T)] 
  print(mito_genes) #可能是13个线粒体基因
  #input_sce=PercentageFeatureSet(input_sce, "^MT-", col.name = "percent_mito")
  input_sce=PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
  fivenum(input_sce@meta.data$percent_mito)
  
  #计算核糖体基因比例
  ribo_genes=rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce),ignore.case = T)]
  print(ribo_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = ribo_genes, col.name = "percent_ribo")
  fivenum(input_sce@meta.data$percent_ribo)
  
  #计算红血细胞基因比例
  Hb_genes=rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce),ignore.case = T)]
  print(Hb_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = Hb_genes,col.name = "percent_hb")
  fivenum(input_sce@meta.data$percent_hb)
  
  #可视化细胞的上述比例情况
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
             "percent_ribo", "percent_hb")
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  p1 
  w=length(unique(input_sce$orig.ident))/3+5;w
  ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
  p2	
  w=length(unique(input_sce$orig.ident))/2+5;w
  ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)
  
  p3=FeatureScatter(input_sce, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
  ggsave(filename="Scatterplot.pdf",plot=p3)
  
  #根据上述指标，过滤低质量细胞/基因
  #过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
  # 一般来说，在CreateSeuratObject的时候已经是进行了这个过滤操作
  # 如果后期看到了自己的单细胞降维聚类分群结果很诡异，就可以回过头来看质量控制环节
  # 先走默认流程即可
  if(F){
    selected_c <- WhichCells(input_sce, expression = nFeature_RNA > 500)
    selected_f <- rownames(input_sce)[Matrix::rowSums(input_sce@assays$RNA$counts > 0 ) > 3]
    input_sce.filt <- subset(input_sce, features = selected_f, cells = selected_c)
    dim(input_sce) 
    dim(input_sce.filt) 
  }
  
  input_sce.filt =  input_sce
  
  # par(mar = c(4, 8, 2, 1))
  # 这里的C 这个矩阵，有一点大，可以考虑随抽样 
  C=subset(input_sce.filt,downsample=100)@assays$RNA$counts
  dim(C)
  C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100

  most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
  
  pdf("TOP50_most_expressed_gene.pdf",width=14)
  boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
          cex = 0.1, las = 1, 
          xlab = "% total count per cell", 
          col = (scales::hue_pal())(50)[50:1], 
          horizontal = TRUE)
  dev.off()
  rm(C)
  
  #过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
  selected_mito <- WhichCells(input_sce.filt, expression = percent_mito < 25)
  selected_ribo <- WhichCells(input_sce.filt, expression = percent_ribo > 3)
  selected_hb <- WhichCells(input_sce.filt, expression = percent_hb < 1 )
  length(selected_hb)
  length(selected_ribo)
  length(selected_mito)
  
  input_sce.filt <- subset(input_sce.filt, cells = selected_mito)
  #input_sce.filt <- subset(input_sce.filt, cells = selected_ribo)
  input_sce.filt <- subset(input_sce.filt, cells = selected_hb)
  dim(input_sce.filt)
  
  table(input_sce.filt$orig.ident) 
  
  #可视化过滤后的情况
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  w=length(unique(input_sce.filt$orig.ident))/3+5;w 
  ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
    NoLegend()
  w=length(unique(input_sce.filt$orig.ident))/2+5;w 
  ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 
  return(input_sce.filt) 
  
}
