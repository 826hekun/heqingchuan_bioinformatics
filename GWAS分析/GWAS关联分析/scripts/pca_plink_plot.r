#!//usr/bin/env Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2020.07.29
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.xet.tech/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.xet.tech/s/26edpc

###############################################################################################

library("argparse")
parser <- ArgumentParser(description='plot PCA scatter from plink')

parser$add_argument( "-i", "--input", type="character",required=T,
		help="input file pc result[required]",
		metavar="filepath")
parser$add_argument("-f", "--groupfile", type="character",required=T,
		help="input group file path[required]",
		metavar="filepath")
parser$add_argument( "-g", "--group", type="character",required=T,
		help="input group name in groupfile file to fill color[required]",
		metavar="group")
parser$add_argument( "-s", "--size", type="integer",required=F,default=3,
		help="point size [optional, default:  %(default)s]",
		metavar="size")
parser$add_argument("-a", "--alpha", type="double", default=0.9,
		help="point transparency [0-1] [optional, default: %(default)s]",
		metavar="alpha")
parser$add_argument("-e", "--ellipse", action='store_true',
		help="whether draw ellipse [optional, default: False]")
parser$add_argument("-L", "--label", action='store_true',
		help="whether show pionts sample name [optional, default:  %(default)s]")
parser$add_argument("-X", "--x.lab", type="character", default="PCA1",
		help="the label for x axis [optional, default:  %(default)s]",
		metavar="label")
parser$add_argument("-Y", "--y.lab", type="character", default="PCA2",
		help="the label for y axis [optional, default:  %(default)s]",
		metavar="label")
parser$add_argument("-T", "--title", type="character", default="PCA",
		help="the label for main title [optional, default:  %(default)s]",
		metavar="label")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default  %(default)s]",
		metavar="path")
parser$add_argument("-n", "--name", type="character", default="demo",
		help="out file name prefix [default  %(default)s]",
		metavar="prefix")
parser$add_argument( "-H", "--height", type="double", default=5,
		help="the height of pic   inches  [default  %(default)s]",
		metavar="number")
parser$add_argument("-W", "--width", type="double", default=6,
		help="the width of pic   inches [default  %(default)s]",
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


library(ggplot2)
library(reshape2)
library(RColorBrewer)


groupfile<-read.table(opt$groupfile,header=T,check.names = F,stringsAsFactors = F,sep="\t",comment.char = "")
pca_table<-read.table(opt$input,header=F,check.names = F,stringsAsFactors = F,fill=T,row.names=1)
#pca_table

#head(groupfile)
#head(pca_table)

rownames(groupfile)<-groupfile[,1]

inter<-intersect(rownames(pca_table),rownames(groupfile))
mdf<-cbind(pca_table[inter,], groupfile[inter,])

mycol=c(brewer.pal(9,"Set1"),brewer.pal(7,"Set2"),brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(7,"Accent"))
shape.value <- c(15:20,1:14,21:25)

head(mdf)
#绘图
p <- ggplot(mdf, aes(x=V3, y=V4, group=get(opt$group)))+
		geom_point(aes(shape=get(opt$group),color=get(opt$group),fill=get(opt$group)),size=opt$size,alpha=opt$alpha) + 
		scale_shape_manual(name=opt$group,values=shape.value)+
		scale_color_manual(name=opt$group,values =mycol)+
		scale_fill_manual(name=opt$group,values =mycol)+
		#stat_ellipse()+
		xlab(opt$x.lab) + ylab(opt$y.lab)+
		labs(title=opt$title)+theme_bw()+ theme(  
				panel.grid=element_blank(), 
				axis.text.x=element_text(colour="black"),
				axis.text.y=element_text(colour="black"),
				panel.border=element_rect(colour = "black"),
				legend.key = element_blank(),plot.title = element_text(hjust = 0.5))
#

# 是否添加置信椭圆
if (opt$ellipse ){
	p = p + stat_ellipse(aes(color=get(opt$group)),level=0.68)
}
if (opt$label){
	p = p + geom_text_repel(label = paste(mdf$V1), colour="black", size=2)
}

pdf(file=paste(opt$outdir,"/",opt$name,".pdf",sep=""), height=opt$height, width=opt$width)
print(p)
dev.off()
png(filename=paste(opt$outdir,"/",opt$name,".png",sep=""), height=opt$height*300, width=opt$width*300, res=300, units="px")
print(p)
dev.off()



