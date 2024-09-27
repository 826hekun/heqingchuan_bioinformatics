#!//usr/bin/env Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2021.10.21
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.xet.tech/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.xet.tech/s/26edpc

###############################################################################################

library("argparse")
parser <- ArgumentParser(description='gwas qq plot')

parser$add_argument( "-i", "--input", type="character",required=T,
		help="input gwas result  file[required]",
		metavar="filepath")
parser$add_argument("-T", "--title", type="character",required=F, default="",
		help="the main title",
		metavar="title")
parser$add_argument( "-c", "--point.color", type="character",required=F,default="blue",
		help="set point color [default=%(default)s]",
		metavar="color")
parser$add_argument("-P", "--point.size", type="double",required=F,default=1,
		help="the  size  point [optional,default:%(default)s]",
		metavar="size")
parser$add_argument("-S", "--point.shape", type="integer",required=F,default=20,
		help="the  shape point ,more info see:https://www.omicsclass.com/article/475 [optional,default:%(default)s]",
		metavar="number")
parser$add_argument("-a", "--axis.size", type="double",required=F,default=7,
		help="the font size of text for axis [optional, default:%(default)s]",
		metavar="size")
parser$add_argument("-r", "--ribbon", action='store_true',
		help="whether draw Confidence interval ribbon    [optional, default: %(default)s]")
parser$add_argument( "-H", "--height", type="double", default=3,
		help="the height of pic   inches  [default %(default)s]",
		metavar="number")
parser$add_argument("-W", "--width", type="double", default=3,
		help="the width of pic   inches [default %(default)s]",
		metavar="number")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default %(default)s]",
		metavar="path")
parser$add_argument("-n", "--name", type="character", default="qq",
		help="out file name prefix [default %(default)s]",
		metavar="prefix")
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}

library(ggplot2)
#qq plot function
gg_qqplot <- function(ps, ci = 0.95,shape=2,size=3,color="black",title="",ribbon=F) {
	n  <- length(ps)
	df <- data.frame(
			observed = -log10(sort(ps)),
			expected = -log10(ppoints(n)),
			clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
			cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
	)
	#log10Pe <- expression(paste("Expected -log"[10], plain(P)))
	log10Pe=expression(Expected~~-log[10](P))
	
	log10Po <- expression(Observed~~-log[10](P))
	p=ggplot(df) +
			geom_point(aes(expected, observed), shape = shape, size = size,color=color) +
			geom_abline(intercept = 0, slope = 1, alpha = 0.5,color="red") +
			# geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
			# geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
			labs(x=log10Pe,y=log10Po,title=title)
	if(ribbon){
		p=p+geom_ribbon(
				mapping = aes(x = expected, ymin = clower, ymax = cupper),
				alpha = 0.1
		)
	}
 	p
}

#read data
data<-read.table(opt$input,header = T,comment.char = "")
colnames(data)<-c("chr","pos","value")

data<-data[!is.nan(data$value),]
#plot qq
p<-gg_qqplot(ps=data$value,shape=opt$point.shape,size=opt$point.size,color=opt$point.color,title=opt$title,ribbon=opt$ribbon) +
		theme_classic()+ theme(  
				#panel.grid=element_blank(), 
				axis.text=element_text(size = opt$axis.size),
				axis.text.x=element_text(colour="black"),
				axis.text.y=element_text(colour="black"),
				#panel.border=element_rect(colour = "black"),
				legend.key = element_blank(),
				plot.title = element_text(hjust = 0.5))+
				scale_y_continuous(expand = c(0,0), limits = c(0, NA))+
				scale_x_continuous(expand = c(0,0), limits = c(0, NA))

png(filename=paste0(opt$outdir,"/",opt$name,".png"), height=opt$height*300, width=opt$width*300, res=300, units="px")
print(p)
dev.off()
pdf(file=paste0(opt$outdir,"/",opt$name,".pdf"), height=opt$height, width=opt$width)
print(p)
dev.off()