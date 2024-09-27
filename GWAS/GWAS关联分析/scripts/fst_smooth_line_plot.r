#!//usr/bin/env Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2021.09.15
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.xet.tech/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.xet.tech/s/26edpc

###############################################################################################

library("argparse")
parser <- ArgumentParser(description='fst smooth line  plot')

parser$add_argument( "-i", "--input", type="character",required=T,
		help="input fst result  file[required]",
		metavar="filepath")
parser$add_argument( "-F", "--fai", type="character",required=T,
		help="genome reference fa index file for filter or order chromosome [required]",
		metavar="fai")
#parser$add_argument("-c", "--threshold", type="double",required=F,default=0.05,
#		help="set top percent threshold as selected region [default 0.05]",
#		metavar="threshold")
parser$add_argument("-s", "--smooth", type="character",required=F,default="DISTANCE",
		help="set smooth and smooth method : DISTANCE or  WINUM   [default %(default)s]",
		metavar="smooth")
parser$add_argument("-N", "--win.num", type="integer",required=F,default=50,
		help="set snp number for   WINUM  smooth method  [default %(default)s]",
		metavar="win.num")
parser$add_argument("-w", "--win.size", type="integer",required=F,default=1000000,
		help="set snp number for   DISTANCE  smooth method   [default %(default)s]",
		metavar="win.size")
parser$add_argument("-f", "--chrlen", type="integer",required=F,default=NULL,
		help="the  length of chromosome threshold to show [default draw all chromosome]",
		metavar="chrlen")

parser$add_argument("-b", "--blank", type="double",required=F,default=0.005,
		help="The gap between the chromosome [optional,default:%(default)s]",
		metavar="blank")

parser$add_argument("-z", "--zscore", action='store_true',
		help="whether zscore fst  [optional, default: False]")
parser$add_argument("-r", "--y.reverse", action='store_true',
		help="whether reverse  Y axis   [optional, default: False]")
parser$add_argument( "-p", "--pallete", type="character",nargs='+',required=F,default=c("#eb65a0", "#22c2e4", "#4abb6b","#f28d21"),
		help="set color pallete  of plot [default=c(\"#eb65a0\", \"#22c2e4\", \"#4abb6b\",\"#f28d21\")]",
		metavar="pallete")

parser$add_argument("-S", "--line.size", type="double",required=F,default=0.5,
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
#parser$add_argument("-Y", "--ylab", type="character", required=F,default="Fst",
#		help="the y axis lab [default %(default)s]",
#		metavar="ylab")
parser$add_argument( "-H", "--height", type="double", default=3,
		help="the height of pic   inches  [default %(default)s]",
		metavar="number")
parser$add_argument("-W", "--width", type="double", default=8,
		help="the width of pic   inches [default %(default)s]",
		metavar="number")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default %(default)s]",
		metavar="path")
parser$add_argument("-n", "--prefix", type="character", default="fst_manhattan",
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

#-----------------------------------------------------------------
# for loess parameter selecting
#-----------------------------------------------------------------

# select span whith min AICc
loess_span_select <- function(value, pos){
	# get spans
	#span <- seq(round(50/length(pos),digits=3), 1, 0.001)
	#span <- seq(round(50/length(pos),digits=2), 1, 0.01)
	#span <- seq(round(50/length(pos),digits=2), 1, 0.03)
	span <- seq(round(50/length(pos),digits=2), 1, 0.05)
	#span <- seq(0.01, 1, 0.05)
	#span <- c(0.025)
	span <- matrix(span,,1)
	cat("span", as.vector(span),"\n")
	# iter
	span.aicc <- apply(span, 1, aicc, value=value, pos=pos)
	cat("span.aicc", as.vector(span.aicc),"\n")
	# get min aicc
	span.notna <- span[!is.na(span.aicc)]
	span.aicc.notna <- span.aicc[!is.na(span.aicc)]
	m <- min(span.aicc.notna)
	if( is.infinite(m) ){
		cat("span selecting failed: min is infinite\n")
		q(status=1)
	}
	# return
	return( span.notna[span.aicc.notna==m] )
}



#-----------------------------------------------------------------
# for smooth method
#-----------------------------------------------------------------
# define function
# for loess
fit_by_loess <- function(x,pos){
	# check the number of points for fitting
	if ( length(x) < 50 ){
		cat("the number of points for fitting is too small(<50)\n")
		q(status=1)
	}
	# get best span
	s <- loess_span_select(x, pos)
	s <- s[1] 	# get one of span
	cat("best_span\t", s, "\n")
	# fitting
	loe <- loess(x ~ pos, span = s, degree=1, family='symmetric', surface='direct')
	# return
	loe$fitted
}
# for window by SNPNUM
win_by_snpnum <- function(x,width){
	p <- x
	len <- length(x)
	for ( i in 1:len ){
		p[i] <- median(x[max(1,i-width):min(len,i+width)])
	}
	p
}
# for window by DISTANCE
# for window by DISTANCE
win_by_distance <- function(x,pos,width){
	p <- x
	len <- length(x)
	for ( i in 1:len ){
		p[i] <- mean(x[ (pos>pos[i]-width) & (pos<pos[i]+width) ])
	}
	p
}
#############################################
#import data
#############################################
rawdata<-read.table(opt$input,header = T,comment.char = "")

rawdata<-rawdata[,c("CHROM", "BIN_START", "BIN_END",    "WEIGHTED_FST")]

colnames(rawdata)<-c("chr","start","end","value")

rawdata$value[rawdata$value<0]=0

#cat("input chromosome name:",unique(rawdata$chr),"\n")

#注意 fai文件中的染色体名称与 GWAS分析中的染色体名称要一致
refl<-read.table(opt$fai,header = F,comment.char = "")
#cat("fa index chr name:",refl$V1,"\n")
if(!is.null(opt$chrlen)){
	refl=subset(refl,V2>=opt$chrlen)
}
rawdata=rawdata[rawdata$chr%in% refl$V1,]


chr=rawdata$chr
value=rawdata$value
pos=rawdata$start
chr_names<-refl$V1
#isremove=rawdata$value
win.value=rawdata$value
## calculate p.value in window
#chr_names <- names( table(chr))
#win.value <- value
#isremove<-value

for ( name in chr_names ) {
	cat("chr_name = ", name, "\n")
	# sub for one chr
	sub <- chr==name
	sub.value <- value[sub]
	sub.pos <- pos[sub]
	# mean in window
	if ( opt$smooth == "WINUM" ){
		win.value[sub] <- win_by_snpnum(x=sub.value, width=opt$win.num)
	} else if ( opt$smooth == "DISTANCE" ){
		win.value[sub] <- win_by_distance(x=sub.value, pos=sub.pos, width=opt$win.size)
	#} else if( opt$smooth == "LOESS" ) {
	#	win.value[sub] <- fit_by_loess(x=sub.value, pos=sub.pos)
	}else{
		stop(paste("not suport smooth method:",opt$smooth))
	}
}


############################
#plot
#############################

#offset chr position
rawdata$smooth.value=win.value
chrlen<-refl$V2
names(chrlen)<-refl$V1
blank=opt$blank*sum(chrlen)
chrpos<-cumsum(chrlen[-length(chrlen)]+blank)
chrpos<-c(0,chrpos)
names(chrpos)<-names(chrlen)
#chrpos<-chrpos-pos[1,names(chrpos)]
chrid<- as.character(rawdata$chr)

ylab="Fst"
quantile(rawdata$value,c(0.01,0.05,0.8,0.90,0.95,0.99),na.rm=T)
threshold=quantile(rawdata$value,1-opt$threshold,na.rm=T)

if(opt$zscore){
	M=mean(rawdata$value)
	S=sd(rawdata$value)
	rawdata$value=(rawdata$value-M)/S
	ylab=expression(Z(Fst))
	#threshold=quantile(rawdata$value,1-opt$threshold)
}
#topdata=subset(rawdata,value>=threshold)
#write.table(topdata,file=paste0(opt$outdir,"/",opt$prefix,"_selected_region_top",opt$threshold,".txt"), quote = F, sep = "\t", row.names = F,)
write.table(rawdata,file=paste0(opt$outdir,"/",opt$prefix,"_","data.txt"), quote = F, sep = "\t", row.names = F,)

#plot manhattan
plot.data<- data.frame(chr=as.character(rawdata$chr),pos=rawdata$start+chrpos[chrid],value=rawdata$smooth.value)
plot.data$chr=factor(plot.data$chr, levels = refl$V1,ordered=T)
p<-ggplot(plot.data,aes(x=pos,y=value,colour=chr))
p<-p+geom_line(size=opt$line.size)
p<-p+scale_colour_manual(values=rep(opt$pallete,times=100))+
		labs(x=opt$xlab,y=ylab,title=opt$title)
xat<- sapply(tapply(plot.data$pos,plot.data$chr,function(x)mean(range(x))), unlist)
p<-p+scale_x_continuous(expand = c(0.001,0),breaks=xat,labels=names(xat)) 
		#geom_hline(yintercept = threshold,linetype = 2, size = 0.2)   #水平阈值线
p<-p+theme_classic()+theme(axis.text=element_text(size = opt$axis.size),
		axis.text.x=element_text(angle = 65,vjust=1,hjust=1,color="black"),  #染色体名称对齐角度
		axis.text.y=element_text(color="black"), 
		axis.title=element_text(size = opt$lab.size), #
		plot.title=element_text(size = opt$title.size,hjust=0.5),
		legend.position='none')+ scale_y_continuous(expand = c(0,0))

if (opt$y.reverse){
	p<-p+scale_y_reverse()
}

png(filename=paste0(opt$outdir,"/",opt$prefix,".png"), height=opt$height*300, width=opt$width*300, res=300, units="px")
print(p)
dev.off()
pdf(file=paste0(opt$outdir,"/",opt$prefix,".pdf"), height=opt$height, width=opt$width)
print(p)
dev.off()


