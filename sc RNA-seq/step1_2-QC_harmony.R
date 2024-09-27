rm(list=ls())
options(stringsAsFactors = F) 
source('scRNA_scripts/lib.R')
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)

###### step1:导入数据 ######   
dir='GSE189357_RAW/' 
fs=list.files('GSE189357_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]

##处理数据，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz到各自的文件夹
#批量将文件名改为 Read10X()函数能够识别的名字
if(F){
lapply(unique(samples),function(x){
  # x = unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE189357_RAW/", paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_"))
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE189357_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE189357_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE189357_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})
}

dir='GSE189357_RAW/'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
}) 
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples  ) 
names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
#看看合并前后的sce变化
sce.all
sce.all <- JoinLayers(sce.all)
sce.all
dim(sce.all[["RNA"]]$counts )

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)
# fivenum(sce.all$nFeature_RNA)
# table(sce.all$nFeature_RNA>800) 
# sce.all=sce.all[,sce.all$nFeature_RNA>800]
# sce.all

library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)

phe$patient = str_split(phe$orig.ident,'[_]',simplify = T)[,2]
table(phe$patient)

###按照患者来源
phe$sample = phe$patient
table(phe$sample)
#肺腺癌（LUAD）从原位腺癌（AIS）到微浸润性腺癌（MIA）和随后的浸润性腺癌（IAC）
phe$sample = gsub("TD5|TD7|TD8", "AIS", phe$sample) 
phe$sample = gsub("TD3|TD4|TD6", "MIA", phe$sample) 
phe$sample = gsub("TD1|TD2|TD9", "IAC", phe$sample) 
table(phe$sample)

sce.all@meta.data = phe


sp='human'
# 如果为了控制代码复杂度和行数 
# 可以省略了质量控制环节
###### step2: QC质控 ######
dir.create("./1-QC")
setwd("./1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
##细胞减少了一点
setwd('../')
getwd()

fivenum(sce.all.filt$percent_ribo)
table(sce.all.filt$nFeature_RNA> 5)

###### step3: harmony整合多个单细胞样品 ######
set.seed(10086)
table(sce.all.filt$orig.ident)
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('../scRNA_scripts/harmony.R')
  # 默认 ScaleData 没有添加"nCount_RNA", "nFeature_RNA"
  # 默认的
  sce.all.int = run_harmony(sce.all.filt)
  setwd('../')
  
}




#######下面代码也可以不运行
###### step4:  看标记基因库 ######
# 原则上分辨率是需要自己肉眼判断，取决于个人经验
# 为了省力，我们直接看 0.1 和 0.8 即可
table(Idents(sce.all.int))
table(sce.all.int$seurat_clusters)
table(sce.all.int$RNA_snn_res.0.1) 
table(sce.all.int$RNA_snn_res.0.2) 
table(sce.all.int$RNA_snn_res.0.8) 

getwd()
dir.create('check-by-0.1')
setwd('check-by-0.1')
sel.clust = "RNA_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 

source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.5')
setwd('check-by-0.5')
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

dir.create('check-by-0.8')
setwd('check-by-0.8')
sel.clust = "RNA_snn_res.0.8"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()

last_markers_to_check
