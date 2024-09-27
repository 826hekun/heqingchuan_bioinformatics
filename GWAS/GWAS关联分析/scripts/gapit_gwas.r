#!//usr/bin/env Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2023.10.07
#version 3.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.xet.tech/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.xet.tech/s/26edpc

###############################################################################################

library("argparse")
parser <- ArgumentParser(description='gwas analysis use gapit')

parser$add_argument( "-i", "--input", type="character",required=T,
		help="input genotype  file , hapmap format [required]",
		metavar="hmp")
parser$add_argument( "-t", "--traits", type="character",required=T,
		help="input traits  file [required] ",
		metavar="traits")
parser$add_argument( "-q", "--structure", type="character",required=F,default=NULL,
		help="input PCA or Q  matrix for structure data [optional ] ",
		metavar="structure")
parser$add_argument( "-k", "--kinship", type="character",required=F,default=NULL,
		help="input kinship  file [optional ] ",
		metavar="kinship")

parser$add_argument( "-m", "--model", type="character",required=F,default="MLM",
		help="model choose for gwas GLM MLM CMLM ECMLM SUPER FarmCPU Blink FaST EMMA EMMAx  [optional default=%(default)s]",
		metavar="model")

parser$add_argument(  "--PCA.total", type="integer",required=F,default=3,
		help="Total Number of PCs as Covariates[optional default=%(default)s]",
		metavar="PCA.total")
#parser$add_argument(  "--PCA.scaling", type="character",required=F,default="None",
#		help="Scale And/Or Center And Scale The SNPs Before Conducting PCA (Scaled, Centered.and.scaled)[optional default=%(default)s]",
#		metavar="PCA.scaling")
#parser$add_argument(  "--kinship.algorithm", type="character",required=F,default="VanRaden",
#		help="Algorithm to Derive Kinship from Genotype (Zhang, Loiselle and EMMA,VanRaden) [optional default=%(default)s]",
#		metavar="kinship.algorithm")
#parser$add_argument(  "--kinship.cluster",nargs='+', type="character",required=F,default="average",
#		help="Clustering algorithm to group individuals based on their kinship (average,complete, ward, single,mcquitty, median, and centroid) [optional default=%(default)s]",
#		metavar="kinship.cluster")
#parser$add_argument(  "--kinship.group",nargs='+', type="character",required=F,default="Mean",
#		help="Method to derive kinship among groups(Mean,Max, Min, and Median) [optional default=%(default)s]",
#		metavar="kinship.group")
#parser$add_argument(  "--group.from", type="integer",required=F,default=1,
#		help="The starting number of groups of iteration[optional default=%(default)s]",
#		metavar="group.from")
#parser$add_argument(  "--group.to", type="integer",required=F,default=1000000,
#		help="The ending number of groups of iteration, if larger than number of individuals (n), n is used[optional default=%(default)s]",
#		metavar="group.to")
#parser$add_argument(  "--group.by", type="integer",required=F,default=10,
#		help="Iteration interval for number of groups[optional default=%(default)s]",
#		metavar="group.by")
parser$add_argument(  "--SNP.fraction", type="double",required=F,default=1,
		help="Fraction of SNPs Sampled to Estimate Kinship and PCs (>0 and <1)[optional default=%(default)s]",
		metavar="SNP.fraction")
parser$add_argument(  "--SNP.MAF", type="double",required=F,default=0,
		help="Minor Allele Frequency to Filter SNPs in GWAS Reports (>0 and <1)[optional default=%(default)s]",
		metavar="SNP.MAF")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [optional default %(default)s]",
		metavar="path")
parser$add_argument("-n", "--name", type="character", default="gwas",
		help="out file name prefix [optional default %(default)s]",
		metavar="prefix")
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}

##########################################################
#GAPTI GWAS analysis
##########################################################


library(multtest)
library("gplots")
#library("LDheatmap")
library("genetics")
library(MASS)
library("compiler")
library(RColorBrewer)
library("scatterplot3d") #Please Download and Install

#source("http://www.zzlab.net/GAPIT/emma.txt")
#source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

if(dir.exists("/share/work/biosoft/GAPIT/")){
  source("/share/work/biosoft/GAPIT/latest/gapit_functions.txt") 
}else if (dir.exists("/biosoft/GAPIT3.0")){

  source("/biosoft/GAPIT3.0/gapit_functions.txt")
  source("/biosoft/GAPIT3.0/emma.txt")
  #source("/biosoft/GAPIT3.0/BLINK.R")
  source("/biosoft/GAPIT3.0/FarmCPU_functions.txt")
  
}else{
  stop("can't find GAPIT source code \n")
}
# if(file.exists("/share/work/biosoft/GAPIT/latest//FarmCPU_functions.txt")){
#   source("/share/work/biosoft/GAPIT/latest//FarmCPU_functions.txt")
#   
# }else{
#   source("/biosoft/GAPIT3.0/FarmCPU_functions.txt")
# }
# if(file.exists("/share/work/biosoft/GAPIT/latest/emma.txt")){
#   source("/share/work/biosoft/GAPIT/latest//emma.txt")
#   
# }else{
#   source("/biosoft/GAPIT3.0/emma.txt")
# }


#source("/biosoft/GAPIT3.0/BLINK.R")




#library(GAPIT3)
#Step 1: Set data directory and import files
myY <- read.table(opt$traits, head = TRUE,stringsAsFactors=F)

myG <- read.table(opt$input, head = FALSE,comment.char="",stringsAsFactors=F)


###get common ID data
samplesID<-intersect(myG[1,],myY[,1])
myY<-myY[myY[,1] %in% samplesID,]
sg<-myG[1,] %in% samplesID
sg[1:11]<-T
myG<-myG[,sg]
#get started

myKI=NULL
myCV=NULL
#myGAPIT=NULL
PCA.total=opt$PCA.total
##Tutorial 3: User defined Kinship and PCs MLM
if(!is.null(opt$kinship) ){
	myKI <- read.table(opt$kinship, head = FALSE,sep="\t",check=F)
	
#	The kinship matrix file (called “KI” in GAPIT) is formatted as an n by n+1 matrix where the first column
#	is the taxa name, and the rest is a square symmetric matrix. Unlike the other input data files, the first row
#	of the kinship matrix file does not consist of headers.
	
#	33-a	2	0.2258837	0.2229
#	22-b	0.2258837	2	0.244
#	22-c	0.2229	0.244	2
}
if(!is.null(opt$structure) ){
	myCV <- read.table(opt$structure, head = FALSE,sep="\t",check=F)
	PCA.total=0
	
#	A file containing covariates (called “CV” in GAPIT) can include information such as population structure
#	(commonly called the “Q matrix”), which are fitted into the GWAS and GS models as fixed effects. These
#	files are formatted similarly to the phenotypic files presented in Section 2.1. Specifically, the first column
#	consists of taxa names, and the remaining columns contain covariate values. The first row consists of
#	column labels. The first column can be labeled “Taxa”, and the remaining columns should be covariate
#	names.
	
#	taxa	Q1	Q2	Q3
#	33-a	0.14	0.972	0.014
#	22-b	0.003	0.993	0.004
	
}






#Methods implimented: 
# 1. GLM (Structure or Q method for GWAS, Pritchard et. al. Genetics, 2000)
# 2. MLM (Q+K, Yu et. al. Nature Genetics, 2006)
# 3. gBLUP (Marker based kinship, Zhang et. al. Journal of Animal Science, 2007)
# 4. PCA (Zhao et. al. Plos Genetics, 2007)
# 5. EMMA (Kang et. al. Genetics, 2008)
# 6. CMLM (Zhang et. al. Nature Genetics, 2010)
# 7. EMMAx (Kang et. al. Nature Genetics, 2010)
# 8. P3D (Zhang et. al. Nature Genetics, 2010)
# 9. FaST-LMM (Lippert et. al. Nature Methods, 2011)
# 10. ECMLM (Li et. al. BMC Bioogy, 2014)
# 11. SUPER (Wang et. al. PLoS One, 2014)

#Literature demonstrated the order of statistical power: BLINK > FarmCPU> MLMM > SUPER >
#		ECMLM > CMLM > MLM > GLM. 

myGAPIT <- GAPIT(
			Y=myY,
			G=myG,
			CV=myCV,
			KI=myKI,
			file.output = FALSE,
			PCA.total=PCA.total,
#			kinship.cluster=opt$kinship.cluster,
#			kinship.group=opt$kinship.group,
#			group.from=opt$group.from,
#			group.to=opt$group.to,
#			group.by=opt$group.by,
			SNP.fraction=opt$SNP.fraction, 
			SNP.MAF=opt$SNP.MAF,
			model=opt$model
	)


write.table(myGAPIT$GWAS,file=paste0(opt$outdir,"/",opt$name,".txt"), quote = F, sep = "\t", row.names = F,)


####################GS ###########################################################
#if(opt$model=="gBLUP"){
#	myGAPIT <- GAPIT(
#			Y=myY,			
#			G=myG,				
#			CV=myCV,
#			KI=myKI,
#			PCA.total=PCA.total,
#			model = "gBLUP",
#			file.output = FALSE
#	#sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
#	#sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
#	#LD=0.1,
#	)
#	write.table(myGAPIT$Pred,file=paste0(opt$outdir,"/",opt$name,".txt"), quote = F, sep = "\t", row.names = F,)
#	
#}
#
#if(opt$model=="cBLUP"){
#	myGAPIT <- GAPIT(
#			Y=myY,			
#			G=myG,				
#			CV=myCV,
#			KI=myKI,
#			PCA.total=PCA.total,
#			model = "cBLUP",
#			file.output = FALSE
#	#sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
#	#sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
#	#LD=0.1,
#	)
#	write.table(myGAPIT$Pred,file=paste0(opt$outdir,"/",opt$name,".txt"), quote = F, sep = "\t", row.names = F,)
#	
#}
#
#
#
#
