# 代码需要保证统一
# 只能是 sce.all


gastric_cancer_markers = c('PTPRC', 
                   'MUC2' , 'ITLN1',
                   'FABP1' , 'APOA1',
                   'CEACAM5' , 'CEACAM6',
                   'EPCAM', 'KRT18', 'MUC1',
                   'MUC6' , 'TFF2',
                   'PGA4' , 'PGA3',
                   'MUC5AC' , 'TFF1','CHGA' , 'CHGB') 
Myo=c("Krt17", "Krt14", "Krt5", "Acta2", "Myl9", "Mylk", "Myh11")
Lum=c("Krt19", "Krt18", "Krt8")
Hs=c("Prlr", "Cited1", "Pgr", "Prom1", "Esr1")  
AV=c("Mfge8", "Trf", "Csn3", "Wfdc18", "Elf5", "Ltf")
Lp=c("Kit", "Aldh1a3", "Cd14")
Fib=c("Col1a1", "Col1a2", "Col3a1", "Fn1")
GSE150580_breast_cancer_markers_list =list(
  Myo=Myo,
  Lum=Lum,
  Hs=Hs, 
  AV=AV,
  Lp=Lp,  
  Fib=Fib 
  
) 

# macrophages (Adgre1, Cd14, and Fcgr3),
# cDCs (Xcr1, Flt3, and Ccr7),
# pDCs (Siglech, Clec10a, and Clec12a), 
# monocytes (Ly6c2 and Spn), 
# neutrophils (Csf3r, S100a8, and Cxcl3),
macrophages=c('Adgre1', 'Cd14',  'Fcgr3')
cDCs=c('Xcr1', 'Flt3',  'Ccr7')
pDCs=c('Siglech', 'Clec10a',  'Clec12a')  
monocytes=c('Ly6c2' , 'Spn')
neutrophils=c('Csf3r', 'S100a8',  'Cxcl3') 
SCP1661_meyloids_markers_list =list(
  macrophages=macrophages,
  cDCs=cDCs,
  pDCs=pDCs, 
  monocytes=monocytes,
  neutrophils=neutrophils  
) 
 
lung_epi_markers =  c('TPPP3',"SPRR3","GDPD3","SPRR1A","SPRR2A","RARRES2","TMPRSS11E",
                    "ASCL3","CFTR","FOXI2","1SG20","FOXI1",
                    "SAA4","SAA2","EFHC1","CCDC153","CCDC113","SAA1","CDC20B","FOXJ1",
                    "MYCL","FOXN4","CCNO",
                    "PIGR","BP1","MUC5A","VMO1","SCGB3A1","CYP2A13","CYP2B6","SCGB1A1",
                    "BCAM","KRT15","KRT5","TP63")

myeloids_markers_list1 =list(
  CM=c("TTN","MYH7","MYH6","TNNT2") ,
  EC=c("VWF", "IFI27", "PECAM1","MGP"),
  FB=c("DCN", "C7" ,"LUM","FBLN1","COL1A2"),
  MP=c("CD163", "CCL4", "CXCL8","PTPRC"),
  SMC=c("ACTA2", "CALD1", "MYH11"),
  Tc=c("CD3D","CD3E"),
  DC1 = c( 'Clec9a', 'Xcr1',   'Wdfy4'), 
  DC2 = c('Itgax', 'Sirpa',   'Cd209a'), 
  mregDCs= c('Ccr7', 'Cd80', 'Cd200',   'Cd247') ,
  hypoxia=c('Hif1a', 'Slc2a1', 'Vegfa', 'Hmox1', 
            'Bnip3', 'Nos2', 'Mmp2', 'Sod3', 
            'Cited2', 'Ldha'),
  peric=c("ABCC9","PDGFRB","RGS5")
)

myeloids_markers_list2 = list(pDC = c("CLEC4C","IRF7","TCF4","GZMB"),
                      cDC1 = c("XCR1","CLNK","CLEC9A"),
                      cDC2 = c("FCER1A","HLA-DPB1","HLA-DQB1","CD1E","CD1C","CLEC10A","HLA-DQA2"),
                      DC3 = c("CCL19","LAMP3","IDO1","IDO2","LAD1","FSCN1","CCR7","LY75","CCL22","CD40","BIRC3","NFKB2"),
                      Macrophages = c("APOC1","HLA-DRB5","C1QA","C1QB"),
                      RTMs = c("THBS1"),#Resident tissue macrophages
                      Lam = c("APOE"),#Lipid associated macrophages
                      Monocytes = c("LYZ","HLA-DRB1","TIMP1","S100A11","CXCL8","IL1B","PTGS2","S100A9","S100A8","MMP19"),
                      Mono_C = c('CD14'),#Mono_CD14
                      Mono_F = c('FCGR3A'),#Mono_FCGR3A
                      Mast = c('TPSAB1' , 'TPSB2'))
 

Tcells_markers = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CCR7', 'SELL' , 'TCF7','CXCR6' , 'ITGA1',
                   'FOXP3', 'IL2RA',  'CTLA4','GZMB', 'GZMK','CCL5',
                   'IFNG', 'CCL4', 'CCL3' ,
                   'PRF1' , 'NKG7') 
###CD4T
CD4_markers_list =list(
  Tc=c("CD3D","CD3E"),
  CD4=c("CD4" ),
  Treg=c("TNFRSF4","BATF","TNFRSF18","FOXP3","IL2RA","IKZF2"),
  naive=c("CCR7","SELL","CD5"),
  Tfh=c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST"),#滤泡辅助性T细胞
  ILC=c("TNFRSF25","KRT81","LST1","AREG","LTB","CD69")
) 

###CD8T
CD8_markers_list1 =list(
  CD8=c("CD8A","CD8B"),
  TN_TCM=c("CCR7","SELL","TCF7","LEF1"),
  TEM=c("GZMK"  ),
  TEFF=c("TBX21","FCGR3A","FGFBP2"),
  TRM=c("XCL1","XCL2","ITGAE","CD69"),
  IEL_T = c("TMIGD2"),
  yT1c=c("GNLY","PTGDS","GZMB","TRDC"),
  yT2c=c("TMN1","HMGB2","TYMS"),
  MAIT_T = c("SLC4A10")
) 
CD8_markers_list2 =list(
  CD8T=c("CD8A","CD8B"),
  MAIT=c("ZBTB16","NCR3","RORA"),
  ExhaustedCD8T=c("LAG3","TIGIT","PDCD1","HAVCR2","CTLA4"),
  EffMemoryCD8=c("EOMES","ITM2C"),
  Resting_NK=c("XCL1","XCL2","KLRC1"),
  Cytotoxic_NK=c("CX3CR1","FGFBP2","FCGR3A","KLRD1"),
  Pre_exhausted=c("IFNG","PRF1","GNLY","GZMA","NKG7","GZMK")
)
 
cd4_and_cd8T_markers_list  =list( 
  naive=c("CCR7","SELL","TCF7","IL7R","CD27","CD28","LEF1","S1PR1"),
  CD8Trm=c("XCL1","XCL2","MYADM"),
  NKTc=c("GNLY","GZMA"), 
  Tfh=c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST"),
  th17=c("IL17A","KLRB1","CCL20","ANKRD28","IL23R","RORC","FURIN","CCR6","CAPG","IL22"),
  CD8Tem=c("CXCR4","GZMH","CD44","GZMK"),
  Treg=c("FOXP3","IL2RA","TNFRSF18","IKZF2"),
  naive=c("CCR7","SELL","TCF7","IL7R","CD27","CD28"),
  CD8Trm=c("XCL1","XCL2","MYADM"), 
  MAIT=c("KLRB1","ZBTB16","NCR3","RORA"),
  yT1c=c("GNLY","PTGDS","GZMB","TRDC"),
  yT2c=c("TMN1","HMGB2","TYMS"),
  yt=c("TRGV9","TRDV2")
) 


# CD20 (MS4A1)表达于除plasma B 之外的所有B，很关键的区分naive 和plasma的marker
# SDC1 = CD138 plasma B （接受抗原，可表达抗体） 
Bcels_markers_list = list(
  All = c('MS4A1','SDC1','CD27','CD38','CD19', 'CD79A'),
  GC_B = c('IL4R','TCL1A','LRMP','SUGCT'),
  IGA_plasm_B= c ( 'IGHA1'), 
  IGG_plasm_B= c ( 'IGHG1')
)  

Hepatic_stellate_markers_list =list(
  qHSC=c("Lrat","Ecm1","Angptl6","Vipr1" ),
  S1=c("Ccl2" ,"Cxcl10" ,"Cxcl1" ,"Ccl7" ),
  S2=c("Acta2" ,"Tpm1" ,"Vim" ,"Tagln","Tnc","Tpm2"),
  S3=c("Col1a1","Col1a2","Col3a1" ,"Lox","Lum" )
)
  
# arteries (HEY1, IGFBP3), capillaries (CD36, CA4), veins (ACKR1) and
# lymphatic ECs (LECs; CCL21, PROX1). 
stromal_markers = c('TEK',"PTPRC","EPCAM","PDPN",
                   "PECAM1",'PDGFRB',"PLVAP",'PROX1','ACKR1','CA4','HEY1',
                   'CSPG4','GJB2', 'RGS5','ITGA7',
                   'ACTA2','RBP1','CD36', 
                   'ADGRE5','COL11A1','FGF7', 'MME') 

last_markers = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'FCGR3A',
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2', ## human  fibo 
                 'GJB2', 'RGS5',
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'MKI67' , 'TOP2A', 
                   'PECAM1', 'VWF',  ## endo 
                 "PLVAP",'PROX1','ACKR1','CA4','HEY1',
                   'EPCAM' , 'KRT19','KRT7', # epi 
                   'FYXD2', 'TM4SF4', 'ANXA4',# cholangiocytes
                   'APOC3', 'FABP1',  'APOA1',  # hepatocytes
                   'Serpina1c',
                   'PROM1', 'ALDH1A1' )

gastric_cancer_markers 
lung_epi_markers
Tcells_markers
stromal_markers 
last_markers 

GSE150580_breast_cancer_markers_list 
SCP1661_meyloids_markers_list 
myeloids_markers_list1 
myeloids_markers_list2 
CD4_markers_list 
CD8_markers_list1 
CD8_markers_list2 
cd4_and_cd8T_markers_list   
Bcels_markers_list 
Hepatic_stellate_markers_list 


markers = c('gastric_cancer_markers','lung_epi_markers',
            'Tcells_markers',
            'stromal_markers', 
            'last_markers' )
markers_list <- c(
  'GSE150580_breast_cancer_markers_list' ,
  'SCP1661_meyloids_markers_list' ,
  'myeloids_markers_list1' ,
  'myeloids_markers_list2' ,
  'CD4_markers_list' ,
  'CD8_markers_list1' ,
  'CD8_markers_list2' ,
  'cd4_and_cd8T_markers_list'   ,
  'Bcels_markers_list' ,
  'Hepatic_stellate_markers_list' 
)

p_umap=DimPlot(sce.all.int, reduction = "umap",raster = F,
               label = T,repel = T) 
p_umap 

if(sp=='human'){
   lapply(markers, function(x){
     #x=markers[1]
     genes_to_check=str_to_upper(get(x)) 
     DotPlot(sce.all.int , features = genes_to_check )  + 
       coord_flip() + 
       theme(axis.text.x=element_text(angle=45,hjust = 1))
     
     h=length( genes_to_check )/6+3;h
     ggsave(paste('check_for_',x,'.pdf'),height = h)
   })
  lapply(markers_list, function(x){
    # x=markers_list[1]
    genes_to_check = lapply(get(x), str_to_upper)
    dup=names(table(unlist(genes_to_check)))[table(unlist(genes_to_check))>1]
    genes_to_check = lapply(genes_to_check, function(x) x[!x %in% dup])
  
    DotPlot(sce.all.int , features = genes_to_check )  + 
     # coord_flip() + 
      theme(axis.text.x=element_text(angle=45,hjust = 1))
    
    w=length( unique(unlist(genes_to_check)) )/5+6;w
    ggsave(paste('check_for_',x,'.pdf'),width  = w)
  })
  
  last_markers_to_check <<- str_to_upper(last_markers ) 

 }else if(sp=='mouse'){
   lapply(markers, function(x){
     #x=markers[1]
     genes_to_check=str_to_title(get(x)) 
     DotPlot(sce.all.int , features = genes_to_check )  + 
       coord_flip() + 
       theme(axis.text.x=element_text(angle=45,hjust = 1))
     
     h=length( genes_to_check )/6+3;h
     ggsave(paste('check_for_',x,'.pdf'),height = h)
   })
   lapply(markers_list, function(x){
     # x=markers_list[1]
     genes_to_check = lapply(get(x), str_to_title)
     dup=names(table(unlist(genes_to_check)))[table(unlist(genes_to_check))>1]
     genes_to_check = lapply(genes_to_check, function(x) x[!x %in% dup])
     
     DotPlot(sce.all.int , features = genes_to_check )  + 
       # coord_flip() + 
       theme(axis.text.x=element_text(angle=45,hjust = 1))
     
     w=length( unique(unlist(genes_to_check)) )/5+6;w
     ggsave(paste('check_for_',x,'.pdf'),width  = w)
   })
   
   last_markers_to_check <<- str_to_title(last_markers ) 
}else {
  print('we only accept human or mouse')
} 

p_all_markers = DotPlot(sce.all.int , features = last_markers_to_check )  + 
  coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1)) 
p_all_markers+p_umap
h=length( last_markers_to_check )/6+3;h
w=length( unique( Idents(sce.all.int)) )/5+10;w
ggsave(paste('last_markers_and_umap.pdf'),width  = w,height = h)

pro = 'qc-'
if("percent_mito" %in% colnames(sce.all.int@meta.data ) ){

  #可视化细胞的上述比例情况
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
  
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1=VlnPlot(sce.all.int , features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  w=length(unique(sce.all.int$orig.ident))/3+5;w
  ggsave(filename=paste0(pro,"Vlnplot1.pdf"),plot=p1,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2=VlnPlot(sce.all.int,  features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
  w=length(unique(sce.all.int$orig.ident))/2+5;w
  ggsave(filename=paste0(pro,"Vlnplot2.pdf"),plot=p2,width = w,height = 5)
  
}
p3=FeatureScatter(sce.all.int , "nCount_RNA", "nFeature_RNA", 
                  pt.size = 0.5)
ggsave(filename=paste0(pro,"Scatterplot.pdf"),plot=p3)


if(T){
  #  remotes::install_github('genecell/COSGR')
  #  genexcell <- Seurat::GetAssayData(object = object[[assay]],slot = slot)
  marker_cosg <- cosg(
    sce.all.int,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100)
  
  save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
  head(marker_cosg)
  
  ## Top10 genes
  library(dplyr)  
  top_10 <- unique(as.character(apply(marker_cosg$names,2,head,10)))
  # width <-0.006*dim(sce.all.int)[2];width
  # height <- 0.25*length(top_10)+4.5;height
  
  width <- 15+0.5*length(unique(Idents(sce.all.int)));width
  height <- 8+0.1*length(top_10);height
  
  sce.Scale <- ScaleData(sce.all.int ,features =  top_10  )  
  
  DoHeatmap(  sce.Scale , top_10 , 
              size=3)
  
  ggsave(filename=paste0(pro,'DoHeatmap_check_top10_markers_by_clusters.pdf') ,
         # limitsize = FALSE,
         units = "cm",width=width,height=height)
  width <- 8+0.6*length(unique(Idents(sce.all.int)));width
  height <- 8+0.2*length(top_10);height
  DotPlot(sce.all.int, features = top_10 ,
          assay='RNA'  )  + coord_flip() +FontSize(y.text = 4)
  ggsave(paste0(pro,'DotPlot_check_top10_markers_by_clusters.pdf'),
         units = "cm",width=width,height=height)
  
  
  ## Top3 genes
  top_3 <- unique(as.character(apply(marker_cosg$names,2,head,3)))
  
  width <- 15+0.2*length(unique(Idents(sce.all.int)));width
  height <- 8+0.1*length(top_3);height
  
  sce.Scale <- ScaleData(sce.all.int ,features =  top_3  )  
  
  DoHeatmap(  sce.Scale , top_3 , 
              size=3)
  ggsave(filename=paste0(pro,'DoHeatmap_check_top3_markers_by_clusters.pdf') ,
         units = "cm",width=width,height=height)
  
  width <- 8+0.2*length(unique(Idents(sce.all.int)));width
  height <- 8+0.1*length(top_3);height
  DotPlot(sce.all.int, features = top_3 ,
          assay='RNA'  )  + coord_flip()
  ggsave(paste0(pro,'DotPlot_check_top3_markers_by_clusters.pdf'),width=width,height=height)
  
}



 
  