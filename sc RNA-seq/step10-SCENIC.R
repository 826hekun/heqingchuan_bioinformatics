#转录因子 (transcription factors, TFs) 是直接作用于基因组，与特定DNA序列结合 (TFBS/motif) ，调控DNA转录过程的一类蛋白质。转录因子可以调节基因组DNA开放性、募集RNA聚合酶进行转录过程、募集辅助因子调节特定的转录阶段，调控诸多生命进程，诸如免疫反应、发育模式等。因此，分析转录因子表达及其调控活性对于解析复杂生命活动具有重要意义。


# 检查并安装所需包
packages <- c("GenomicFeatures","AUCell", "RcisTarget", "GENIE3", "zoo", "mixtools", "rbokeh", 
              "DT", "NMF", "pheatmap", "R2HTML", "Rtsne", "doMC", "doRNG")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# 对于从github安装的包，使用devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("SCopeLoomR", quietly = TRUE)) {
  devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
}


library(GenomicFeatures)
library(AUCell)
library(RcisTarget)
library(GENIE3)

# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)
library(R2HTML)
library(Rtsne)
library(doMC)
library(doRNG)
library(usethis)


if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")


#人类
#load in the motifannotation this will load it into your environment but as the name in which is given to the list argument
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")

#rename the motif annnotion by attributing it to the variable that is in the error
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

set.seed(123)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)


##SCENIC（single-cell regulatory network inference and clustering）是一个基于共表达和motif分析，计算单细胞转录组数据基因调控网络重建以及细胞状态鉴定的方法。在输入单细胞基因表达量矩阵后，SCENIC经过以下三个步骤完成转录因子分析：

#第一步，GENIE3（随机森林)/GRNBoost (Gradient Boosting) 推断转录因子与候选靶基因之间的共表达模块，每个模块包含一个转录因子及其靶基因，纯粹基于共表达；
#第二步，RcisTatget分析每个共表达模块中的基因，以鉴定enriched motifs，仅保留TF motif富集的模块和targets，构建TF-targets网络，每个TF及其潜在的直接targets gene被称作一个调节因子（Regulons）；
#第三步，AUCelll计算调节因子（Regulons）的活性，这将确定Regulon在哪些细胞中处于“打开”状态。


##==分析准备==##
dir.create("10-SCENIC")
setwd("10-SCENIC/") 
dir.create("int")
scRNA=readRDS( "../7-epi/epi_sce_celltype.rds")
table(scRNA$celltype)
cancer = scRNA[, scRNA$celltype %in% c( 'Cancer(CRABP2+)','Cancer(TM4SF1+)', 'Cancer(UBE2C+)','Clara-like cancer')]

##准备细胞meta信息
cellInfo <- data.frame(cancer@meta.data)
# colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
# colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
# #colnames(cellInfo)[which(colnames(cellInfo)=="celltype_Monaco")] <- "celltype"

cellInfo <- cellInfo[,c("sample","celltype")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
#为了节省计算资源，随机抽取1000个细胞的数据子集
# subcell <- sample(colnames(cancer),300)
# scRNAsub <- cancer[,subcell]
# table( scRNAsub$celltype )
Idents(cancer) = cancer$celltype
allCells=names(Idents(cancer))
allType = levels(Idents(cancer))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(cancer)== x ]
  cg=sample(cgCells,50)
  cg
}))
scRNAsub = cancer[, allCells %in% choose_Cells]
saveRDS(scRNAsub, "scRNAsub.rds")

exprMat <- as.matrix(scRNAsub@assays$RNA$counts)
dim(exprMat)

##设置分析环境
getwd()
#根据对应物种下载相应的数据库文件到目录下。可以使用多种下载方式，这里用的是下载完成后本地上传。
#下载链接https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/
mydbDIR <- "/home/data/t050453/paper_plot_redraw/scRNA/GSE189357-LUAD-scRNA-ST/10-SCENIC"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
           "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp", "10kb")
scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=60,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "SCENIC")
saveRDS(scenicOptions, "int/scenicOptions.rds")

##==转录调控网络推断==##
##基因过滤
#过滤标准是基因表达量之和>细胞数*3%，且在1%的细胞中表达
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##计算相关性矩阵
runCorrelation(exprMat_filtered, scenicOptions)
##TF-Targets相关性回归分析
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
#这一步消耗的计算资源非常大，这里300个基因运行不到20分钟，还行


##推断共表达模块
#TF是转录因子名称，Target是潜在靶基因的名字，weight是TF与Target之间的相关性权重。
runSCENIC_1_coexNetwork2modules(scenicOptions)

##推断转录调控网络（regulon）
runSCENIC_2_createRegulons(scenicOptions)
##这一步运行很慢
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量

scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=1,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "erythrogenesis")
##==regulon活性评分与可视化==##
##regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(scRNA@assays$RNA$counts)
exprMat_all <- log2(exprMat_all+1)
saveRDS(exprMat_all,"exprMat_all.rds")
# rm(list=ls())
# scenicOptions<-readRDS("scenicOptions.rds")
# exprMat_all<-readRDS("exprMat_all.rds")
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

#runSCENIC_3_scoreCells(scenicOptions, exprMat=log2(as.matrix(scRNA@assays$RNA@counts)+1))



# #使用shiny互动调整阈值
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
# savedSelections <- shiny::runApp(aucellApp)
# #保存调整后的阈值
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)
savehistory("scenic_code.txt")


regulons <- loadInt(scenicOptions, "regulons")

#一个regulon可以对应多个靶标，而AUC则是在靶标大于10个的regulons中进行运算。
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))


regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
head(tableSubset)


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
##随机选取20个
a = as.data.frame(regulonActivity_byCellType)
a <- a[sample(rownames(a),20),]
regulonActivity_byCellType_Scaled <- t(scale(t(a), center = T, scale=T))

dev.new()
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()

##细胞特异性调控因子
#我们还可以计算每个regulon在细胞中的特异性系数（Regulon Specificity Score, RSS）来衡量regulon在不同细胞类型间的特异程度。
# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
cellAnnotation=cellInfo[colnames(regulonAUC), "celltype"]
cellAnnotation = na.omit(cellAnnotation)
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation)
rssPlot <- plotRSS(rss)
print(rssPlot$plot)
