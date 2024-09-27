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
parser <- ArgumentParser(description='rod manhattan plot')

parser$add_argument( "-i", "--input", type="character",required=T,
		help="input rod result  file[required]",
		metavar="filepath")
parser$add_argument( "-F", "--fai", type="character",required=T,
		help="genome reference fa index file for filter or order chromosome [required]",
		metavar="fai")
parser$add_argument("-c", "--threshold", type="double",required=F,default=0.05,
		help="set top percent threshold as selected region [default %(default)s]",
		metavar="threshold")
parser$add_argument("-f", "--chrlen", type="integer",required=F,default=NULL,
		help="the  length of chromosome threshold to show [default draw all chromosome]",
		metavar="chrlen")

parser$add_argument("-b", "--blank", type="double",required=F,default=0.005,
		help="The gap between the chromosome [optional,default:%(default)s]",
		metavar="blank")

parser$add_argument("-z", "--zscore", action='store_true',
		help="whether zscore ROD  [optional, default: %(default)s]")
parser$add_argument("-r", "--y.reverse", action='store_true',
		help="whether reverse  Y axis   [optional, default: %(default)s]")
parser$add_argument("-v", "--vline", action='store_true',
		help="plot vline instead of point [optional, default: %(default)s]")
parser$add_argument( "-p", "--pallete", type="character",nargs='+',required=F,default=c("#eb65a0", "#22c2e4", "#4abb6b","#f28d21"),
		help="set color pallete  of plot [default=c(\"#eb65a0\", \"#22c2e4\", \"#4abb6b\",\"#f28d21\")]",
		metavar="pallete")
parser$add_argument("-P", "--point.size", type="double",required=F,default=1,
		help="the point size [optional,default:%(default)s]",
		metavar="size")
parser$add_argument("-S", "--line.size", type="double",required=F,default=0.2,
		help="the line  width [optional,default:%(default)s]",
		metavar="size")

parser$add_argument("-l", "--lab.size", type="double",required=F,default=8,
		help="the font size of x and y label [optional, default:%(default)s]",
		metavar="size")
parser$add_argument("-t", "--title.size", type="double",required=F,default=12,
		help="the title size [optional, default: %(default)s]",
		metavar="size")
parser$add_argument("-a", "--axis.size", type="double",required=F,default=7,
		help="the font size of text for axis [optional, default:%(default)s]",
		metavar="size")

parser$add_argument("-T", "--title", type="character",required=F, default="",
		help="the main title",
		metavar="title")

parser$add_argument("-X", "--xlab", type="character", required=F,default="Chromosome",
		help="the x axis lab [default %(default)s]",
		metavar="label")
parser$add_argument("-Y", "--ylab", type="character", required=F,default="ROD",
		help="the y axis lab [default %(default)s]",
		metavar="label")
parser$add_argument( "-H", "--height", type="double", default=3,
		help="the height of pic   inches  [default %(default)s]",
		metavar="number")
parser$add_argument("-W", "--width", type="double", default=8,
		help="the width of pic   inches [default %(default)s]",
		metavar="number")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default %(default)s]",
		metavar="path")
parser$add_argument("-n", "--prefix", type="character", default="rod_manhattan",
		help="out file name prefix [default %(default)s]",
		metavar="prefix")
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}

require(ggplot2)
col=c("#4197d8", "#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d")
#import data
rawdata<-read.table(opt$input,header = T,comment.char = "")

rawdata<-rawdata[,c("CHROM", "BIN_START", "BIN_END",    "ROD")]

colnames(rawdata)<-c("chr","start","end","value")

#cat("input chromosome name:",unique(rawdata$chr),"\n")

#注意 fai文件中的染色体名称与 染色体名称要一致
refl<-read.table(opt$fai,header = F,comment.char = "")
#cat("fa index chr name:",refl$V1,"\n")
if(!is.null(opt$chrlen)){
	refl=subset(refl,V2>=opt$chrlen)
}
rawdata=rawdata[rawdata$chr%in% refl$V1,]

#offset  chr position
chrlen<-refl$V2
names(chrlen)<-refl$V1
blank=opt$blank*sum(chrlen)
chrpos<-cumsum(chrlen[-length(chrlen)]+blank)
chrpos<-c(0,chrpos)
names(chrpos)<-names(chrlen)
#chrpos<-chrpos-pos[1,names(chrpos)]
chrid<- as.character(rawdata$chr)

ylab="ROD"
quantile(rawdata$value,c(0.01,0.05,0.8,0.90,0.95,0.99),na.rm=T)
threshold=quantile(rawdata$value,1-opt$threshold,na.rm=T)

if(opt$zscore){
	M=mean(rawdata$value)
	S=sd(rawdata$value)
	rawdata$value=(rawdata$value-M)/S
	ylab=expression(Z(ROD))
	threshold=quantile(rawdata$value,1-opt$threshold)
}
topdata=subset(rawdata,value>=threshold)
write.table(topdata,file=paste0(opt$outdir,"/",opt$prefix,"_selected_region_top",opt$threshold,".txt"), quote = F, sep = "\t", row.names = F,)

#plot manhattan
plot.data<- data.frame(chr=as.character(rawdata$chr),pos=rawdata$start+chrpos[chrid],value=rawdata$value)

#去掉小于0的区域
plot.data=subset(plot.data,value>=0)

plot.data$chr=factor(plot.data$chr, levels = refl$V1,ordered=T)

p<-ggplot(plot.data,aes(x=pos,y=value,colour=chr))

if (opt$vline){
	p<-p+  geom_segment(aes(x=pos, xend=pos, y=0, yend=value,colour=chr), 
			size=opt$line.size) 
}else{
	p<-p+geom_point(size=opt$point.size)
}

p<-p+scale_colour_manual(values=rep(opt$pallete,times=100))+
		labs(x=opt$xlab,y=ylab,title=opt$title)

xat<- sapply(tapply(plot.data$pos,plot.data$chr,function(x)mean(range(x))), unlist)

p<-p+scale_x_continuous(expand = c(0.001,0),breaks=xat,labels=names(xat)) 
		if(opt$vline){
			#p<-p+geom_hline(yintercept = threshold,linetype = 2, size = 0.2)   #水平阈值线
		}else{
			p<-p+geom_hline(yintercept = threshold,linetype = 2, size = 0.2)   #水平阈值线
		}
p<-p+theme_classic()+theme(axis.text=element_text(size = opt$axis.size),
		axis.text.x=element_text(angle = 65,vjust=1,hjust=1,color="black"),  #染色体名称对齐角度
		axis.text.y=element_text(color="black"), 
		axis.title=element_text(size = opt$lab.size), #
		plot.title=element_text(size = opt$title.size,hjust=0.5),
		legend.position='none')+ scale_y_continuous(expand = c(0,0))

if (opt$y.reverse){
	p<-p+scale_y_reverse()
}

png(filename=paste0(opt$outdir,"/",opt$prefix,"_manhattan.png"), height=opt$height*300, width=opt$width*300, res=300, units="px")
print(p)
dev.off()
pdf(file=paste0(opt$outdir,"/",opt$prefix,"_manhattan.pdf"), height=opt$height, width=opt$width)
print(p)
dev.off()
