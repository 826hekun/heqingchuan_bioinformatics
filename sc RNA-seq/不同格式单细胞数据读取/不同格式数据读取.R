#清空环境 加载需要的R包
rm(list=ls())
options(stringsAsFactors = F) 
source('./lib.R')

##10X标准格式
dir='GSE212975_10x/'
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
sce.all <- JoinLayers(sce.all)

dim(sce.all[["RNA"]]$counts )

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 



##h5格式
#清空环境 加载需要的R包
rm(list=ls())
options(stringsAsFactors = F) 
source('./lib.R')

library(hdf5r)
library(stringr)
library(data.table)

dir='GSE215120_h5/'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X_h5(file.path(dir,pro )) 
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
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 




#RDS文件
rm(list=ls())
options(stringsAsFactors = F) 
source('scRNA_scripts/lib.R')

dir='outputs/' 
samples=list.files( dir,pattern = 'rds',full.names = T,recursive = T )
samples 

library(data.table)
sceList = lapply(samples,function(pro){ 
  #pro=samples[1] 
  print(pro) 
  ct=readRDS(pro)  
  ct[1:4,1:4]
  sce=CreateSeuratObject(  ct , 
                           project =  gsub('_gex_raw_counts.rds','',
                                           basename(pro) ) ,
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
sce.all <- JoinLayers(sce.all)


##txt.gz格式
#清空环境 加载需要的R包
rm(list=ls())
options(stringsAsFactors = F) 
source('./lib.R')


dir='GSE167297_txt/'
samples=list.files( dir )
samples 

sceList = lapply(samples,function(pro){ 
    # pro=samples[1] 
    print(pro)  
    ct=fread(file.path( dir ,pro),data.table = F)
    ct[1:4,1:4]
    rownames(ct)=ct[,1]
    colnames(ct) = paste(gsub('_CountMatrix.txt.gz','',pro),
                         colnames(ct) ,sep = '_')
    ct=ct[,-1] 
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
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 

##csv.gz格式
#清空环境 加载需要的R包
rm(list=ls())
options(stringsAsFactors = F) 
source('./lib.R')


dir='GSE129516_csv/'
samples=list.files( dir )
samples 

sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  colnames(ct) = paste(gsub('_filtered_gene_bc_matrices.csv.gz','',pro),
                       colnames(ct) ,sep = '_')
  ct=ct[,-1] 
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
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 

#读取三个独立的文件
#加载需要的R包
rm(list=ls())
options(stringsAsFactors = F) 
source('./lib.R')

library(data.table)
library(Matrix)
# https://hbctraining.github.io/scRNA-seq/lessons/readMM_loadData.html
#本次的样品为smart-seq2的样品
#将三个文件按照对应的格式分别读取进来
mtx=readMM( "./GSE184708/GSE184708_raw_counts_gonad_all_samples_XX_XY_E10_to_E16.mtx.gz" )
mtx[1:4,1:4]
dim(mtx)

cl=fread( "./GSE184708/GSE184708_mayere_barcodes.tsv.gz" ,
          header = F,data.table = F )
head(cl)

rl=fread( "./GSE184708/GSE184708_mayere_genes.tsv.gz" ,
          header = F,data.table = F )
head(rl) 

#整合矩阵信息
colnames(mtx)=cl$V1
rownames(mtx)=rl$V1

#创建seurat对象
sce.all=CreateSeuratObject(counts = mtx ,
                           project = 'mouse',
                           min.cells = 5,
                           min.features = 300 )

#进行分组
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident)
library(stringr)
phe=str_split(rownames(sce.all@meta.data),'_',simplify = T)
head(phe)
table(phe[,2])
sce.all$group=phe[,2]
table(phe[,3])
sce.all$sex=phe[,3]
table(phe[,1])
sce.all$orig.ident=phe[,1]