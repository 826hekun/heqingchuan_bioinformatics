#!/biosoft/miniconda/bin/Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2020.07.29
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.h5.xeknow.com/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.h5.xeknow.com/s/26edpc

###############################################################################################


library("argparse")
parser <- ArgumentParser(description='plot tree nwk file')

parser$add_argument( "-t", "--tree", type="character",required=T,
		help="phylip tree file nwk format [required]",
		metavar="filepath")

 parser$add_argument( "-l", "--layout", type="character",required=F,default="circular",
 		help="set tree layout:unrooted,circular,slanted,rectangular [default %(default)s]",
 		metavar="layout")
parser$add_argument(  "--outgroup", type="character",required=F,default=NULL,nargs="+",
                     help="set tree outgroup [default %(default)s]",
                     metavar="outgroup")

parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default %(default)s]",
		metavar="path")
parser$add_argument("-p", "--prefix", type="character", default="tree",
		help="out file name prefix [default %(default)s]",
		metavar="prefix")
parser$add_argument( "-H", "--height", type="double", default=10,
		help="the height of pic   inches  [default %(default)s]",
		metavar="number")
parser$add_argument("-W", "--width", type="double", default=10,
		help="the width of pic   inches [default %(default)s]",
		metavar="number")

opt <- parser$parse_args()

#parser$print_help()

#set some reasonable defaults for the options that are needed,
#but were not specified.
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}





library(ggtree)
library(ggplot2)
library(ape)
library("RColorBrewer")
clist <- list(
		"shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
		"strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
		"oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
		"keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
		"vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
		"muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
		"teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
		"merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
		"funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
		"retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
		"cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
		"cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
		"morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
		"wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
		"krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))


tree <- read.tree(opt$tree)

if(!is.null(opt$outgroup)){
  
  tree <- root(tree, which(tree$tip.label %in% opt$outgroup))
}

#opt$layout="rectangular"
    p=ggtree(tree, ladderize=FALSE, size=0.3, branch.length="none",layout=opt$layout)+ geom_tiplab(size=4)+
      geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 1))
      
  
  
pdf(file=paste(opt$outdir, "/",opt$prefix,".pdf", sep=""),h=opt$height,w=opt$width)
print(p)
dev.off()
png(file=paste(opt$outdir, "/",opt$prefix,".png" ,sep=""),h=opt$height*300,w=opt$width*300,res=300)
print(p)
dev.off()

