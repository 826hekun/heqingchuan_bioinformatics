#!//usr/bin/env Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2023.10.07
#version 1.1
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.xet.tech/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.xet.tech/s/26edpc

###############################################################################################

library("argparse")
parser <- ArgumentParser(description='fst and pi scattor plot')

parser$add_argument( "-F", "--fst", type="character",required=T,
		help="input pop1 vs pop2  fst result  file[required]",
		metavar="FstFile")

parser$add_argument( "-a", "--pi1", type="character",required=T,
		help="input pop1 pi result  file[required]",
		metavar="Pi1File")

parser$add_argument( "-b", "--pi2", type="character",required=T,
		help="input pop2 pi result  file[required]",
		metavar="Pi2File")

parser$add_argument( "-A", "--pop1", type="character",required=T,
		help="set pop1 group name result  file[required]",
		metavar="pop1name")

parser$add_argument( "-B", "--pop2", type="character",required=T,
		help="set pop2 group name[required]",
		metavar="pop2name")

parser$add_argument("-c", "--threshold", type="double",required=F,default=0.05,
		help="set top percent threshold as selected region [default %(default)s]",
		metavar="threshold")

parser$add_argument("-z", "--zscore", action='store_true',
		help="whether  using zscore transformation of  fst value [optional, default: False]")
parser$add_argument("-L", "--log2", action='store_true',
		help="whether using log2 transformation of   pi ratio  [optional, default: False]")

parser$add_argument( "-p", "--pallete", type="character",nargs='+',required=F,default=c("#4197d8", "#f8c120", "#413496", "#495226"),
		help="set color pallete  of plot [default=c(\"#4197d8\", \"#f8c120\", \"#413496\", \"#495226\")]",
		metavar="color")
parser$add_argument("-P", "--point.size", type="double",required=F,default=0.2,
		help="the point size [optional,default:%(default)s]",
		metavar="size")
parser$add_argument("-S", "--line.size", type="double",required=F,default=0.2,
		help="the line  width [optional,default:%(default)s]",
		metavar="size")

parser$add_argument("-l", "--lab.size", type="double",required=F,default=8,
		help="the font size of x and y label [optional, default:%(default)s]",
		metavar="size")
parser$add_argument("-t", "--title.size", type="double",required=F,default=12,
		help="the point size [optional, default: %(default)s]",
		metavar="size")
parser$add_argument("-x", "--axis.size", type="double",required=F,default=7,
		help="the font size of text for axis [optional, default:%(default)s]",
		metavar="size")

parser$add_argument("-T", "--title", type="character",required=F, default="",
		help="the main title",
		metavar="title")
parser$add_argument("-f", "--format", type="character",required=F, default="pdf",
		help="the output pic format, png or pdf [default %(default)s]",
		metavar="format")
parser$add_argument("-X", "--xlab", type="character", required=F,default="Chromosome",
		help="the x axis lab [default %(default)s]",
		metavar="label")
parser$add_argument("-Y", "--ylab", type="character", required=F,default="Fst",
		help="the y axis lab [default %(default)s]",
		metavar="label")
parser$add_argument( "-H", "--height", type="double", default=8,
		help="the height of pic   inches  [default %(default)s]",
		metavar="number")
parser$add_argument("-W", "--width", type="double", default=8,
		help="the width of pic   inches [default %(default)s]",
		metavar="number")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default %(default)s]",
		metavar="path")
parser$add_argument("-n", "--name", type="character", default="fst_and_pi",
		help="out file name prefix [default %(default)s]",
		metavar="prefix")
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}




fst<-read.table(opt$fst,sep="\t",header=T,comment.char="")

fst<-data.frame(CHROM=fst$CHROM,BIN_START=fst$BIN_START, BIN_END=fst$BIN_END,FST=fst$WEIGHTED_FST)

#fst$ID=paste(fst$CHROM,fst$BIN_START, fst$BIN_END)
pop1<-read.table(opt$pi1,sep="\t",header=T,comment.char="")
pop1<-data.frame(CHROM=pop1$CHROM,BIN_START=pop1$BIN_START, BIN_END=pop1$BIN_END,PIA=pop1$PI)

#pop1$ID=paste(pop1$CHROM,pop1$BIN_START, pop1$BIN_END)
pop2<-read.table(opt$pi2,sep="\t",header=T,comment.char="")
pop2<-data.frame(CHROM=pop2$CHROM,BIN_START=pop2$BIN_START, BIN_END=pop2$BIN_END,PIB=pop2$PI)

#pop2$ID=paste(pop2$CHROM,pop2$BIN_START, pop2$BIN_END)

d<-merge(fst,pop1,by=c( "CHROM" ,"BIN_START",  "BIN_END"))
data<-merge(d,pop2,by=c( "CHROM" ,"BIN_START",  "BIN_END"))

data$PIR=data$PIA/data$PIB



#xlab=expression(log[2](pi[ bquote(.(labNames[1]) ]/pi[",opt$pop2,"]"))
#ylab=expression(paste(Z(Fst), "$popA vs $popB"))

xlab= bquote(pi[.(opt$pop1)]/pi[.(opt$pop2)])
ylab= expression(F[st])



fst_cut=0
pi_up_cut=0
pi_down_cut=0

if(opt$zscore){
	M=mean(data$FST)
	S=sd(data$FST)
	data$FST=(data$FST-M)/S
	ylab=expression(z(F[st]))
	cat("z Fst quantile:\n",quantile(data$FST,c(0.8,0.90,0.95,0.99)),"\n")
	
	fst_cut=quantile(data$FST,1-opt$threshold)
	fst_legend=c(bquote(z(F[st])<.(round(fst_cut,2))),
			bquote(z(F[st])>=.(round(fst_cut,2))),
			expression("Cumulative (%)"))
	
	
}else{
	cat("Fst quantile:\n",quantile(data$FST,c(0.8,0.90,0.95,0.99)),"\n")
	
	fst_cut=quantile(data$FST,1-opt$threshold)
	
	
	fst_legend=c(bquote(F[st]<.(round(fst_cut,2))),
			bquote(F[st]>=.(round(fst_cut,2))),
			expression("Cumulative (%)"))
}

if(opt$log2){
	xlab=bquote(log[2](pi[.(opt$pop1)]/pi[.(opt$pop2)]))
	data$PIR=log2(data$PIR)
	cat("log2 Pi ratio(popa/popb) quantile:\n",quantile(data$PIR,c(0.8,0.90,0.95,0.99)),"\n")
	
	pi_up_cut=quantile(data$PIR,1-opt$threshold)
	pi_down_cut=quantile(data$PIR,opt$threshold)
	pi_legend=c(substitute(paste(down,"<", log2(pi[ratio]),"<",up),
					list(up=round(pi_up_cut,2),down=round(pi_down_cut,2))),
			bquote(log2(pi[ratio])<=.(round(pi_down_cut,2))),
			bquote(log2(pi[ratio])>=.(round(pi_down_cut,2))),
			expression("Cumulative (%)"))			
}else{
	cat("Pi ratio(popa/popb) quantile:\n",quantile(data$PIR,c(0.8,0.90,0.95,0.99)),"\n")
	
	pi_up_cut=quantile(data$PIR,1-opt$threshold)
	pi_down_cut=quantile(data$PIR,opt$threshold)
	
	pi_legend=c(substitute(paste(down,"<", pi[ratio],"<",up),
					list(up=round(pi_up_cut,2),down=round(pi_down_cut,2))),
			bquote(pi[ratio]<=.(round(pi_down_cut,2))),
			bquote(pi[ratio]>=.(round(pi_up_cut,2))),
			expression("Cumulative (%)"))
	
	
}

#set color 
library(RColorBrewer)
mycol<-brewer.pal(8,"Set1")
palette(c(mycol[3],"#000000BB",mycol[2],"grey",mycol[5]))


#palette(c("black","grey","blue",mycol[4]))
data$is_select=NA
data$is_select[data$FST>=fst_cut & data$PIR>=pi_up_cut]=1    #pop B  受选择区域
data$is_select[data$FST>=fst_cut & data$PIR<=pi_down_cut]=3  # pop A  受选择区域
data$is_select[is.na(data$is_select)]=2
#fst_cut
#pi_up_cut
#pi_down_cut

#输出受选择区域文件：1：pop B  受选择区域   3：pop A  受选择区域
#注意 此文件如果做了转换，即为转换后的数据

out.df<-subset(data,is_select==1 | is_select==3)
out.df$is_select[out.df$is_select==1]=paste0(opt$pop2," selected region")
out.df$is_select[out.df$is_select==3]=paste0(opt$pop1," selected region")
out.df$is_select[out.df$is_select==2]=paste0("not selected region")

write.table(out.df,
		file=paste0(opt$outdir,"/",opt$name,".selected.region.txt"), 
		sep = "\t",row.names = F,quote =F)

write.table(data,
		file=paste0(opt$outdir,"/",opt$name,".merged.all.txt"), 
		sep = "\t",row.names = F,quote =F)

if (opt$format=="pdf"){
	#dev.off()
	pdf(file=paste0(opt$outdir,"/",opt$name,".pdf"),w=opt$width,h=opt$height)
}else{
	#dev.off()
	png(filename=paste0(opt$outdir,"/",opt$name,".png"),w=opt$width*300,h=opt$height*300,res=300)
}

#分割画布
layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), widths=c(5,2),heights=c(2,5), TRUE)
#Make Scatterplot
par(mar=c(5.1,4.1,0,0))
plot(data$FST~data$PIR,pch=21,
		col=data$is_select,bg=data$is_select,
		xlim=c(floor(min(data$PIR)),ceiling(max(data$PIR))),
		ylim=c(floor(min(data$FST)),ceiling(max(data$FST))),
		xlab=xlab,
		ylab=ylab,
		cex=0.5)
legend("bottomright",legend=c(paste0(" Selected region (",opt$pop1,")"),
				paste0(" Selected region (",opt$pop2,")")),
		pch=21,col=c(3,1), pt.bg=c(3,1),bty="n",pt.cex=1.5, adj = c(0, 0.5))
abline(h=fst_cut,col=2,lty=2,lwd=1.5)
abline(v=pi_up_cut,col=2,lty=2,lwd=1.5)
abline(v=pi_down_cut,col=2,lty=2,lwd=1.5)


#up pi bar 
par(mar=c(0,4.1,3,0))
xhist <- hist(data$PIR,breaks=seq(from=floor(min(data$PIR)),to=ceiling(max(data$PIR)),length.out=100),
		plot=F)
#设置pi柱子颜色向量

xbar.col=rep(4,100)
xbar.col[xhist$mids<pi_down_cut]=3
xbar.col[xhist$mids>pi_up_cut]=1

xhist <- hist(data$PIR,breaks=seq(from=floor(min(data$PIR)),to=ceiling(max(data$PIR)),length.out=100),
		ann=FALSE,axes=T,xaxt="n",col=xbar.col,border=xbar.col)
#up line Cumulative
ec <- ecdf(data$PIR)
lines(x = xhist$mids, y=ec(xhist$mids)*max(xhist$counts), col ='black',lwd=2)
axis(4, at=seq(from = 0, to = max(xhist$counts), length.out = 6), labels=seq(0, 100, 20), col = 'black', col.axis = 'black')
mtext(side = 4, line = 3, 'Cumulative  (%)', col = 'black',cex=0.8,adj = 1)
mtext(side = 2, line = 3, 'Counts', col = 'black',cex=0.8)
abline(v=pi_up_cut,col=2,lty=2,lwd=1.5)
abline(v=pi_down_cut,col=2,lty=2,lwd=1.5)
legend("bottomright",legend=pi_legend,
		lty=1,lwd=1.5,col=c(4,3,1,"black"),bty="n")

#right bar
ymin=floor(min(data$FST))
ymax=ceiling(max(data$FST))
range=ymax-ymin

par(mar=c(5.1,0,0,3))
yhist <- hist(data$FST,breaks=seq(from=floor(min(data$FST)),to=ceiling(max(data$FST)),length.out=100),
		plot=F)

#bar颜色向量设置
#设置柱子颜色向量

ybar.col=rep(4,100)
ybar.col[yhist$mids>fst_cut]=5

#水平bar绘图
ec <- ecdf(data$FST)
barplot(yhist$counts,horiz=TRUE,space=0,axes=T,col=ybar.col,main="",border=ybar.col)

#right line Cumulative
axis(3, at=seq(from = 0, to = max(yhist$counts), length.out = 6), labels=seq(0, 100, 20), col = 'black', col.axis = 'black')
lines(y = (yhist$mids-ymin)/range*100, x=ec(yhist$mids)*max(yhist$counts), col ='black',lwd=2)
mtext(side = 3, line = 3, 'Cumulative (%)', col = 'black',cex=0.8,padj=1,adj = 1)
mtext(side = 1, line = 3, 'Counts', col = 'black',cex=0.8)
abline(h=(fst_cut-ymin)/range*100-0.5,col=2,lty=2,lwd=1.5)

legend("bottomright",legend=fst_legend,
		lty=1,lwd=1.5,col=c(4,5,"black"),bty="n")

dev.off()



